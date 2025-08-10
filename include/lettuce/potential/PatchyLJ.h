#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

// template<typename T>
// const std::vector<T> patchAngles = {0, 2*M_PI/4, 4*M_PI/4, 6*M_PI/4};
// const std::vector<T> patchAngles = {0, 2*M_PI/6, 2*M_PI/6, 3*M_PI/6, 4*M_PI/6, 5*M_PI/6};

template<typename T>
inline T wrapAngle(const T& ang) {
    return std::fmod(ang + M_PI, static_cast<T>(2.0) * M_PI) - M_PI; 
}

template<typename T>
T findBestPatchSignedAngle(const T& phi, const T& gamma, const std::vector<T>& patchAngs) {
    T minAbs = M_PI;
    T bestSigned = 0;
    for (unsigned int k = 0; k < patchAngs.size(); ++k) {
        T patchAngle = phi + patchAngs[k];
        T deltaTheta = wrapAngle(patchAngle - gamma);
        T absDelta = std::abs(deltaTheta);
        if (absDelta < minAbs) {
            minAbs = absDelta;
            bestSigned = deltaTheta;
        }
    }
    return bestSigned;
}

template<typename T>
inline T findPatchAngleUnwrapped(const T& phi, const T& gamma, const T& patchAngle) {
    return phi + patchAngle - gamma;
}

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJPotential {
public:
    PatchyLJPotential(T epsilon, T sigma, T cutoff, T sigmaAngular, std::vector<T> patchAngles)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngular(sigmaAngular),
          m_sigma6(std::pow(sigma, 6)),
          m_2sigmaAng_sq(2.0 * sigmaAngular * sigmaAngular),
          m_patchAngles(std::move(patchAngles)) {}

T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
    if (R <= static_cast<T>(0) || R > m_cutoff) return static_cast<T>(0);

    const T LJpot = lennardJones(R, m_sigma, m_epsilon);

    if (R < m_sigma) {
        return LJpot;
    }

    const auto dx = dr[0];
    const auto dy = dr[1];

    const T phi_i = p1.phi;
    const T phi_j = p2.phi;

    const auto& patchAngs = m_patchAngles;
    const std::size_t patchNum_i = patchAngs.size();
    const std::size_t patchNum_j = patchNum_i;

    const T gamma_ij = std::atan2(dy, dx);
    const T gamma_ji = wrapAngle(gamma_ij + M_PI);

    if (patchNum_i == 0 || patchNum_j == 0) return static_cast<T>(0);

    T sum_ei = static_cast<T>(0);
    T sum_ej = static_cast<T>(0);
    for (std::size_t l = 0; l < patchNum_i; ++l) {
        const T theta_i = wrapAngle(findPatchAngleUnwrapped(phi_i, gamma_ij, patchAngs[l]));
        sum_ei += std::exp(-(theta_i * theta_i) / m_2sigmaAng_sq);
    }
    for (std::size_t m = 0; m < patchNum_j; ++m) {
        const T theta_j = wrapAngle(findPatchAngleUnwrapped(phi_j, gamma_ji, patchAngs[m]));
        sum_ej += std::exp(-(theta_j * theta_j) / m_2sigmaAng_sq);
    }

    const T A = sum_ei * sum_ej;

    return LJpot * A;
}

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_sigma6, m_2sigmaAng_sq;
    std::vector<T> m_patchAngles;
};



template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJForce {
public:
    PatchyLJForce(T epsilon, T sigma, T cutoff, T sigmaAngular, std::vector<T> patchAngles)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngular(sigmaAngular),
          m_sigma6(std::pow(sigma, 6)),
          m_2sigmaAng_sq(static_cast<T>(2.0) * sigmaAngular * sigmaAngular),
          m_patchAngles(std::move(patchAngles)) {}

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= static_cast<T>(0) || R > m_cutoff) return {{static_cast<T>(0), static_cast<T>(0), static_cast<T>(0)}};

        const T dx = dr[0];
        const T dy = dr[1];
        const T r2 = dx * dx + dy * dy;

        const T LJpot = lennardJones(R, m_sigma, m_epsilon);
        const T dVdR = lennardJonesDerivative(R, m_sigma, m_epsilon);

        const T dRdx = dx / R;
        const T dRdy = dy / R;
        const T dVdx = dVdR * dRdx;
        const T dVdy = dVdR * dRdy;

        const T fx_lj = -dVdx;
        const T fy_lj = -dVdy;

        if (R < m_sigma) {
            return {{fx_lj, fy_lj, 0}};
        }

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;

        const auto& patchAngs = m_patchAngles;
        const std::size_t patchNum_i = patchAngs.size();
        const std::size_t patchNum_j = patchNum_i;

        const T gamma_ij = std::atan2(dy, dx);
        const T gamma_ji = wrapAngle(gamma_ij + M_PI);

        if (patchNum_i == 0 || patchNum_j == 0) {
            return {{fx_lj, fy_lj, static_cast<T>(0)}};
        }

        T sum_ei = static_cast<T>(0);
        T sum_ej = static_cast<T>(0);
        T Si = static_cast<T>(0);
        T Sj = static_cast<T>(0);

        for (std::size_t l = 0; l < patchNum_i; ++l) {
            const T theta_i = wrapAngle(findPatchAngleUnwrapped(phi_i, gamma_ij, patchAngs[l]));
            const T ei = std::exp(-(theta_i * theta_i) / m_2sigmaAng_sq);
            sum_ei += ei;
            Si += (static_cast<T>(-2.0) * theta_i / m_2sigmaAng_sq) * ei;
        }
        for (std::size_t m = 0; m < patchNum_j; ++m) {
            const T theta_j = wrapAngle(findPatchAngleUnwrapped(phi_j, gamma_ji, patchAngs[m]));
            const T ej = std::exp(-(theta_j * theta_j) / m_2sigmaAng_sq);
            sum_ej += ej;
            Sj += (static_cast<T>(-2.0) * theta_j / m_2sigmaAng_sq) * ej;
        }

        // A and its theta-derivative sums
        const T A = sum_ei * sum_ej;
        // total derivative of A wrt theta_i summed over all pairs: sum_ej * Si
        const T dA_dtheta_i_total = sum_ej * Si;
        // total derivative of A wrt theta_j summed over all pairs: sum_ei * Sj
        const T dA_dtheta_j_total = sum_ei * Sj;

        // derivatives of gamma_ij wrt x_j (particle j coordinates)
        const T dgammadx = -dy / r2;   // ∂γ_ij/∂x_j
        const T dgammady = dx / r2;    // ∂γ_ij/∂y_j

        // ∂A/∂x_j = ∂A/∂θ_i * ∂θ_i/∂γ_ij * ∂γ_ij/∂x_j
        //         + ∂A/∂θ_j * ∂θ_j/∂γ_ji * ∂γ_ji/∂γ_ij * ∂γ_ij/∂x_j
        // Here ∂θ_i/∂γ_ij = -1, ∂θ_j/∂γ_ji = -1, and ∂γ_ji/∂γ_ij = 1
        const T dA_dx = (dA_dtheta_i_total * (static_cast<T>(-1.0))
                          + dA_dtheta_j_total * (static_cast<T>(-1.0))) * dgammadx;
        const T dA_dy = (dA_dtheta_i_total * (static_cast<T>(-1.0))
                          + dA_dtheta_j_total * (static_cast<T>(-1.0))) * dgammady;

        // Total force = - (dV/dx * A + LJpot * dA/dx) etc.
        const T fx = -(dVdx * A + LJpot * dA_dx) ;
        const T fy = -(dVdy * A + LJpot * dA_dy) ;

        // Torque on j: -∂U/∂φ_j = -LJpot * dA/dθ_j_total
        const T torquej = -LJpot * dA_dtheta_j_total;

        // add LJ-only force components (note: above already includes LJ part)
        // We computed fx,fy including LJ contribution; if you'd rather split,
        // return fx_lj + angular part etc. Current formula matches original's
        // combined expression per pair.

        return {{fx, fy, torquej}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_sigma6, m_2sigmaAng_sq;
    std::vector<T> m_patchAngles;
};