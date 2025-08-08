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
constexpr T wrapAngle(const T& ang) {
    return std::fmod(ang + M_PI, 2.0 * M_PI) - M_PI; 
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
T findPatchAngle(const T& phi, const T& gamma, const T& patchAngle) {
    return wrapAngle(phi + patchAngle - gamma);
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
    if (R <= 0 || R > m_cutoff) return 0;

    const T LJpot = lennardJones(R, m_sigma, m_epsilon);

    if (R < m_sigma) {
        return LJpot;
    }

    const auto dx = dr[0];
    const auto dy = dr[1];

    const T phi_i = p1.phi;
    const T phi_j = p2.phi;

    const std::vector<T> phi_i_patchAngs = m_patchAngles;
    const std::vector<T> phi_j_patchAngs = m_patchAngles;
    const T patchNum_i = phi_i_patchAngs.size();
    const T patchNum_j = phi_j_patchAngs.size();

    const T gamma_ij = std::atan2(dy, dx);
    const T gamma_ji = wrapAngle(gamma_ij + M_PI);

    T A = 0;

    for (int l = 0; l < patchNum_i; l++){
        for (int m = 0; m < patchNum_j; m++){
            T theta_i = findPatchAngle(phi_i, gamma_ij, phi_i_patchAngs[l]);
            T theta_j = findPatchAngle(phi_j, gamma_ji, phi_j_patchAngs[m]);
            A += std::exp(-theta_i * theta_i / m_2sigmaAng_sq) *
                 std::exp(-theta_j * theta_j / m_2sigmaAng_sq);
        }
    }

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
          m_2sigmaAng_sq(2.0 * sigmaAngular * sigmaAngular),
          m_patchAngles(std::move(patchAngles)) {}

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};

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

        const std::vector<T> phi_i_patchAngs = m_patchAngles;
        const std::vector<T> phi_j_patchAngs = m_patchAngles;
        const T patchNum_i = phi_i_patchAngs.size();
        const T patchNum_j = phi_j_patchAngs.size();

        const T gamma_ij = std::atan2(dy, dx);
        const T gamma_ji = wrapAngle(gamma_ij + M_PI);

        T A, fx, fy, torquej = 0;

        for (int l = 0; l < patchNum_i; l++){
            for (int m = 0; m < patchNum_j; m++){

                // θ_i = φ_i + φ_k - γ_ij
                // θ_j = φ_j + φ_k - (γ_ji + M_PI)
                T theta_i = findPatchAngle(phi_i, gamma_ij, phi_i_patchAngs[l]);
                T theta_j = findPatchAngle(phi_j, gamma_ji, phi_j_patchAngs[m]);

                A = std::exp(-theta_i * theta_i / m_2sigmaAng_sq) *
                            std::exp(-theta_j * theta_j / m_2sigmaAng_sq);

                // ============== FORCE CALCULATION ==============
                T dgammadx = -dy / r2;   // ∂γ_ij/∂x_j
                T dgammady = dx / r2;    // ∂γ_ij/∂y_j

                // Derivatives of A w.r.t. gamma_ij
                T dA_dtheta_i = A * (-2.0 * theta_i / m_2sigmaAng_sq);
                T dA_dtheta_j = A * (-2.0 * theta_j / m_2sigmaAng_sq);
                
                // ∂A/∂x_j = ∂A/∂θ_i * ∂θ_i/∂γ_ij * ∂γ_ij/∂x_j 
                //         + ∂A/∂θ_j * ∂θ_j/∂γ_ji * ∂γ_ji/∂γ_ij * ∂γ_ij/∂x_j
                T dA_dx = dA_dtheta_i * (-1.0) * dgammadx 
                                + dA_dtheta_j * (-1.0) * 1.0 * dgammadx;
                T dA_dy = dA_dtheta_i * (-1.0) * dgammady 
                                + dA_dtheta_j * (-1.0) * 1.0 * dgammady;

                fx += -(dVdx * A + LJpot * dA_dx);
                fy += -(dVdy * A + LJpot * dA_dy);

                // ============== TORQUE CALCULATION ==============
                torquej += -LJpot * dA_dtheta_j;  // -∂U/∂φ_j = -LJpot * ∂A/∂θ_j
            }
        }
        return {{fx, fy, torquej}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_sigma6, m_2sigmaAng_sq;
    std::vector<T> m_patchAngles;
};