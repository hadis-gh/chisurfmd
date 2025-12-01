#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

std::vector<std::vector<double>> patchMat;

template<typename T>
inline T wrapAngle(const T& ang) {
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
inline T findPatchAngleUnwrapped(const T& phi, const T& gamma, const T& patchAngle) {
    return phi + patchAngle - gamma;
}

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJPotential {
public:
    PatchyLJPotential(T epsilon, T sigma, T cutoff, T sigmaAngular, int patchNums)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngular(sigmaAngular),
          m_sigma6(std::pow(sigma, 6)),
          m_2sigmaAng_sq(2.0 * sigmaAngular * sigmaAngular),
          m_patchNums(patchNums)
    {
        // Generate evenly spaced patch angles
        m_patchAngles.reserve(m_patchNums);
        for (int a = 0; a < m_patchNums; a++) {
            m_patchAngles.push_back(2 * M_PI * a / m_patchNums); 
        }
        
        // Initialize patchMat to all ones if not provided or incorrect size
        // if (patchMat.empty() || patchMat.size() != static_cast<std::size_t>(m_patchNums)) {
        //     patchMat.resize(m_patchNums);
        //     for (int i = 0; i < m_patchNums; i++) {
        //         patchMat[i].resize(m_patchNums, 1.0);
        //     }
        // }
    }

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return 0;

        const T baseLJ = lennardJones(R, m_sigma, m_epsilon);

        if (R < m_sigma) {
            return baseLJ;
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

        if (patchNum_i == 0 || patchNum_j == 0) return 0;

        T sum_ei = 0;
        T sum_ej = 0;
        for (std::size_t l = 0; l < patchNum_i; ++l) {
            const T theta_i = wrapAngle(findPatchAngleUnwrapped(phi_i, gamma_ij, patchAngs[l]));
            sum_ei += std::exp(-(theta_i * theta_i) / m_2sigmaAng_sq);
        }
        for (std::size_t m = 0; m < patchNum_j; ++m) {
            const T theta_j = wrapAngle(findPatchAngleUnwrapped(phi_j, gamma_ji, patchAngs[m]));
            sum_ej += std::exp(-(theta_j * theta_j) / m_2sigmaAng_sq);
        }

        const T A = sum_ei * sum_ej;

        return baseLJ * A;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_sigma6, m_2sigmaAng_sq;
    int m_patchNums;
    std::vector<T> m_patchAngles;
};



template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJForce {
public:
    PatchyLJForce(T epsilon, T sigma, T cutoff, T sigmaAngular, int patchNums)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngular(sigmaAngular),
          m_sigma6(std::pow(sigma, 6)),
          m_2sigmaAng_sq(2.0 * sigmaAngular * sigmaAngular),
          m_patchNums(patchNums) 
    {
        // Generate evenly spaced patch angles
        m_patchAngles.reserve(m_patchNums);
        for (int a = 0; a < m_patchNums; a++) {
            m_patchAngles.push_back(2 * M_PI * a / m_patchNums); 
        }
        // Initialize patchMat to all ones if not provided or incorrect size
        // if (patchMat.empty() || patchMat.size() != static_cast<std::size_t>(m_patchNums)) {
        //     patchMat.resize(m_patchNums);
        //     for (int i = 0; i < m_patchNums; i++) {
        //         patchMat[i].resize(m_patchNums, 1.0);
        //     }
        // }
    }

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};

        const T dx = dr[0];
        const T dy = dr[1];
        const T r2 = dx * dx + dy * dy;

        const T baseLJ = lennardJones(R, m_sigma, m_epsilon);
        
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
            return {{fx_lj, fy_lj, 0}};
        }

        T sum_ei = 0;
        T sum_ej = 0;
        T Si = 0;
        T Sj = 0;

        for (std::size_t l = 0; l < patchNum_i; ++l) {
            const T theta_i = wrapAngle(findPatchAngleUnwrapped(phi_i, gamma_ij, patchAngs[l]));
            const T ei = std::exp(-(theta_i * theta_i) / m_2sigmaAng_sq);
            sum_ei += ei;
            Si += (-2.0 * theta_i / m_2sigmaAng_sq) * ei;
        }
        for (std::size_t m = 0; m < patchNum_j; ++m) {
            const T theta_j = wrapAngle(findPatchAngleUnwrapped(phi_j, gamma_ji, patchAngs[m]));
            const T ej = std::exp(-(theta_j * theta_j) / m_2sigmaAng_sq);
            sum_ej += ej;
            Sj += (-2.0 * theta_j / m_2sigmaAng_sq) * ej;
        }

        const T A = sum_ei * sum_ej;
        const T dA_dtheta_i_total = sum_ej * Si;
        const T dA_dtheta_j_total = sum_ei * Sj;

        const T dgammadx = -dy / r2;   // ∂γ_ij/∂x_j
        const T dgammady = dx / r2;    // ∂γ_ij/∂y_j

        const T dA_dx = -(dA_dtheta_i_total + dA_dtheta_j_total) * dgammadx;
        const T dA_dy = -(dA_dtheta_i_total + dA_dtheta_j_total) * dgammady;

        const T fx = -(dVdx * A + baseLJ * dA_dx);
        const T fy = -(dVdy * A + baseLJ * dA_dy);
        const T torque_j = -baseLJ * dA_dtheta_j_total;
        const T torque_i = -baseLJ * dA_dtheta_i_total;

        return {{fx, fy, torque_j}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_sigma6, m_2sigmaAng_sq;
    int m_patchNums;
    std::vector<T> m_patchAngles;
};