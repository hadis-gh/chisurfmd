#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

template<typename T>
constexpr T wrapAngle(const T& ang) {
    return std::fmod(ang + M_PI, 2.0 * M_PI) - M_PI; 
}

template<typename T>
T findBestPatchAngle(const T& phi, const T& gamma, const std::vector<T>& patchAngs) {
    T minTheta = M_PI;
    for (unsigned int k = 0; k < patchAngs.size(); ++k) {
        T patchAngle = phi + patchAngs[k];
        T deltaTheta = wrapAngle(patchAngle - gamma);
        minTheta = std::min(minTheta, std::abs(deltaTheta));
    }
    return minTheta;
}
template<typename T>
T findBestPatchAngle(const T& phi, const T& gamma, const std::vector<T>& patchAngs, T& bestPatchAngle, T& theta) {
    T minTheta = M_PI;
    for (auto patchAng : patchAngs) {
        T angle = wrapAngle(phi + patchAng - gamma);
        T absAngle = std::abs(angle);
        if (absAngle < minTheta) {
            minTheta = absAngle;
            bestPatchAngle = phi + patchAng;
            theta = angle;
        }
    }
    return minTheta;
}

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJPotential {
public:
    PatchyLJPotential(T epsilon, T sigma, T cutoff, T sigmaAngular)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngular(sigmaAngular),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)),
          m_2sigmaAng_sq(2.0 * sigmaAngular * sigmaAngular) {}

T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
    if (R <= 0 || R > m_cutoff) return 0;

    const T sr6 = m_sigma6 / std::pow(R, 6);
    const T sr12 = sr6 * sr6;
    const T LJpot = 4.0 * m_epsilon * (sr12 - sr6);

    if (R < m_sigma) {
        return LJpot;
    }

    const T phi_i = p1.phi;
    const T phi_j = p2.phi;
    // const std::vector<T> phi_i_patcheAngs = p1.patchAngs;
    // const std::vector<T> phi_j_patcheAngs = p2.patchAngs;

    const std::vector<T> phi_i_patcheAngs {0, M_PI/2, M_PI, 3* M_PI/2};
    const std::vector<T> phi_j_patcheAngs {0, M_PI/2, M_PI, 3* M_PI/2};

    const T gamma_ij = std::atan2(dr[1], dr[0]);
    const T gamma_ji = wrapAngle(gamma_ij + M_PI);

    T minTheta_i = findBestPatchAngle(phi_i, gamma_ij, phi_i_patcheAngs);
    T minTheta_j = findBestPatchAngle(phi_j, gamma_ji, phi_j_patcheAngs);

    const T angularFactor = std::exp(-minTheta_i * minTheta_i / m_2sigmaAng_sq) *
                             std::exp(-minTheta_j * minTheta_j / m_2sigmaAng_sq);

    return LJpot * angularFactor;
}

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_sigma6, m_sigma12, m_2sigmaAng_sq;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJForce {
public:
    PatchyLJForce(T epsilon, T sigma, T cutoff, T sigmaAngular)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngular(sigmaAngular),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)),
          m_2sigmaAng_sq(2.0 * sigmaAngular * sigmaAngular) {}

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;

        const std::vector<T> patchAngs_i {0, M_PI/2, M_PI, 3*M_PI/2};
        const std::vector<T> patchAngs_j {0, M_PI/2, M_PI, 3*M_PI/2};

        // Lennard-Jones terms
        const T sr6 = m_sigma6 / std::pow(R, 6);
        const T sr12 = sr6 * sr6;
        const T LJ_pot = 4.0 * m_epsilon * (sr12 - sr6);
        const T dLJ_dR = 24.0 * m_epsilon / R * (2.0 * sr12 - sr6);

        if (R < m_sigma) {
            const T force_mag = dLJ_dR;
            return {{force_mag * dr[0] / R, force_mag * dr[1] / R, 0}};
        }

        // Angles
        const T gamma_ij = std::atan2(dr[1], dr[0]);
        const T gamma_ji = wrapAngle(gamma_ij + M_PI);

        T bestPatchAngle_i, theta_kij;
        T theta_i = findBestPatchAngle(phi_i, gamma_ij, patchAngs_i, bestPatchAngle_i, theta_kij);

        T bestPatchAngle_j, theta_kji;
        T theta_j = findBestPatchAngle(phi_j, gamma_ji, patchAngs_j, bestPatchAngle_j, theta_kji);

        const T angularFactor = std::exp(-theta_i * theta_i / m_2sigmaAng_sq) *
                                std::exp(-theta_j * theta_j / m_2sigmaAng_sq);

        // Derivatives of angular factor
        const T dA_dtheta_i = -theta_i / m_2sigmaAng_sq * angularFactor;
        const T dA_dtheta_j = -theta_j / m_2sigmaAng_sq * angularFactor;

        // ∂theta/∂gamma = -1 ⇒ ∂theta/∂x = -∂gamma/∂x
        const T dx = dr[0], dy = dr[1];
        const T R2 = R * R;

        const T dgamma_dxi = -dy / R2;
        const T dgamma_dyi = dx / R2;

        const T dtheta_i_dx = -dgamma_dxi;
        const T dtheta_i_dy = -dgamma_dyi;
        const T dtheta_j_dx = dgamma_dxi;
        const T dtheta_j_dy = dgamma_dyi;

        const T dA_dx = dA_dtheta_i * dtheta_i_dx + dA_dtheta_j * dtheta_j_dx;
        const T dA_dy = dA_dtheta_i * dtheta_i_dy + dA_dtheta_j * dtheta_j_dy;

        // Derivative w.r.t. phi_i (only θ_i depends on it)
        const T dtheta_i_dphi_i = 1.0;
        const T dA_dphi_i = dA_dtheta_i * dtheta_i_dphi_i;

        // Final force and torque
        const T fx = -(dLJ_dR * angularFactor * dx / R + LJ_pot * dA_dx);
        const T fy = -(dLJ_dR * angularFactor * dy / R + LJ_pot * dA_dy);
        const T torque = -LJ_pot * dA_dphi_i;

        return {{fx, fy, torque}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_sigma6, m_sigma12, m_2sigmaAng_sq;
};