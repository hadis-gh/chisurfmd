#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

template<typename T>
const std::vector<T> patchAngles = {0};
// const std::vector<T> patchAngles = {0, M_PI/3, 2*M_PI/3, 3*M_PI/3, 4*M_PI/3, 5*M_PI/3};

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

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJPotential {
public:
    PatchyLJPotential(T epsilon, T sigma, T cutoff, T sigmaAngular)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngular(sigmaAngular),
          m_sigma6(std::pow(sigma, 6)),
          m_2sigmaAng_sq(2.0 * sigmaAngular * sigmaAngular) {}

T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
    if (R <= 0 || R > m_cutoff) return 0;

    const T sr6 = m_sigma6 / std::pow(R, 6);
    const T sr12 = sr6 * sr6;
    const T LJpot = 4.0 * m_epsilon * (sr12 - sr6);

    if (R < m_sigma) {
        return LJpot;
    }

    const auto dx = dr[0];
    const auto dy = dr[1];

    const T phi_i = p1.phi;
    const T phi_j = p2.phi;

    const std::vector<T> phi_i_patcheAngs = patchAngles<T>;
    const std::vector<T> phi_j_patcheAngs = patchAngles<T>;

    const T gamma_ij = std::atan2(dy, dx);
    const T gamma_ji = wrapAngle(gamma_ij + M_PI);

    T minTheta_i = findBestPatchAngle(phi_i, gamma_ij, phi_i_patcheAngs);
    T minTheta_j = findBestPatchAngle(phi_j, gamma_ji, phi_j_patcheAngs);

    const T angularFactor = std::exp(-minTheta_i * minTheta_i / m_2sigmaAng_sq) *
                            std::exp(-minTheta_j * minTheta_j / m_2sigmaAng_sq);

    return LJpot * angularFactor;
}

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_sigma6, m_2sigmaAng_sq;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJForce {
public:
    PatchyLJForce(T epsilon, T sigma, T cutoff, T sigmaAngular)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngular(sigmaAngular),
          m_sigma6(std::pow(sigma, 6)),
          m_2sigmaAng_sq(2.0 * sigmaAngular * sigmaAngular) {}

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};

        const auto dx = dr[0];
        const auto dy = dr[1];
        const auto r2 = dx * dx + dy * dy;
        const auto r6 = r2 * r2 * r2;

        const T sr6 = m_sigma6 / std::pow(R, 6);
        const T sr12 = sr6 * sr6;
        const T LJpot = 4.0 * m_epsilon * (sr12 - sr6);

        const T fx_lj = - 4.0 * m_epsilon 
            * (6.0 * m_sigma6 / (r6 * R) - 12.0 * m_sigma6 * m_sigma6 / (r6 * r6 * R)) * dx / R; 

        const T fy_lj = - 4.0 * m_epsilon 
            * (6.0 * m_sigma6 / (r6 * R) - 12.0 * m_sigma6 * m_sigma6 / (r6 * r6 * R)) * dy / R; 

        if (R < m_sigma) {
            return {{fx_lj, fy_lj, 0}};
        }

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;

        const std::vector<T> phi_i_patcheAngs = patchAngles<T>;
        const std::vector<T> phi_j_patcheAngs = patchAngles<T>;

        const T gamma_ij = std::atan2(dy, dx);
        const T gamma_ji = wrapAngle(gamma_ij + M_PI);

        T minTheta_i = findBestPatchAngle(phi_i, gamma_ij, phi_i_patcheAngs);       // θi = φi + patch - γ
        T minTheta_j = findBestPatchAngle(phi_j, gamma_ji, phi_j_patcheAngs);       // θj = φj + patch - (γ + π)
        
        const T angularFactor = std::exp(-minTheta_i * minTheta_i / m_2sigmaAng_sq) *
                        std::exp(-minTheta_j * minTheta_j / m_2sigmaAng_sq);

        const T dAdphi = - 2.0 * minTheta_j * angularFactor / m_2sigmaAng_sq;

        const T fx = fx_lj * angularFactor - dAdphi * dy / r2 * LJpot;          // -∂U/∂xj = - A * ∂U/∂x - ∂A/∂x * U = fx_lj - ∂A/∂θ * ∂θ/∂γ * ∂γ/dx * Ulj

        const T fy = fy_lj * angularFactor + dAdphi * dx / r2 * LJpot;

        const T torquei = - dAdphi * LJpot;                                     // -∂U/∂φj = -∂A/∂φj * Ulj

        return {{fx, fy, torquei}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_sigma6, m_2sigmaAng_sq;
};