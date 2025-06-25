#pragma once

#include <vector>
#include <cmath>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

template<typename T>
constexpr T wrapAngle(const T& ang){
    return std::fmod(ang + M_PI, 2*M_PI) - M_PI; 
}

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJPotential {
public:
    PatchyLJPotential(T epsilon, T sigma, T cutoff, T sigmaAngularScale)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngularScale(sigmaAngularScale),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)) {}

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R == 0 || R > m_cutoff) return 0;

        T LJpot = 4 * m_epsilon * (
            m_sigma12 / std::pow(R, 12) - 
            m_sigma6 / std::pow(R, 6)
        );

        if (R <= m_epsilon) {
            return LJpot;
        }

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T gamma = std::atan2(dr[1], dr[0]);
        const int patchNum = 3;

        T angularTerm = 0;

        for (int n = 0; n < patchNum; ++n) {
            const T patchAngle = n* 2 * M_PI /patchNum;

            const T tetha_i = wrapAngle(phi_i + patchAngle - gamma);
            const T tetha_j = wrapAngle(phi_j + patchAngle - gamma);

            angularTerm += std::exp(-tetha_i*tetha_i/(2*m_sigmaAngularScale)) * 
                           std::exp(-tetha_j*tetha_j/(2*m_sigmaAngularScale));
        }

        return LJpot * angularTerm;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngularScale;
    T m_sigma6, m_sigma12;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJForce {
public:
    PatchyLJForce(T epsilon, T sigma, T cutoff, T sigmaAngularScale)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_sigmaAngularScale(sigmaAngularScale),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)) {}

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R == 0 || R > m_cutoff) return {{0, 0}};

        // Lennard-Jones potential and its derivative
        T LJpot = 4 * m_epsilon * (
            m_sigma12 / std::pow(R, 12) - 
            m_sigma6 / std::pow(R, 6)
        );
        T dVLJ_dR = 4 * m_epsilon * (
            -12 * m_sigma12 / std::pow(R, 13) + 6 * m_sigma6 / std::pow(R, 7)
        );

        // Angular part and its derivatives
        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T gamma = std::atan2(dr[1], dr[0]);
        const int patchNum = 3;

        T angularTerm = 0;
        T dAngular_dPhi = 0;
        T dAngular_dx = 0;
        T dAngular_dy = 0;

        for (int n = 0; n < patchNum; ++n) {
            const T patchAngle = n* 2 * M_PI /patchNum;
            const T tetha_i = wrapAngle(phi_i + patchAngle - gamma);
            const T tetha_j = wrapAngle(phi_j + patchAngle - gamma);
            const T exp_i = std::exp(-tetha_i*tetha_i/(2*m_sigmaAngularScale));
            const T exp_j = std::exp(-tetha_j*tetha_j/(2*m_sigmaAngularScale));
            const T term = exp_i * exp_j;
            const auto drr2 = dr[0]*dr[0] + dr[1]*dr[1];

            angularTerm += term;
            dAngular_dPhi += -2* tetha_i* term / (2 * m_sigmaAngularScale);
            dAngular_dx += (-2* tetha_i* exp_i *(-dr[1]/drr2)*exp_j + -2* tetha_j* exp_j *(-dr[1]/drr2)*exp_i) / (2 * m_sigmaAngularScale) ;
            dAngular_dy += (-2* tetha_i* exp_i *(dr[0]/drr2)*exp_j + -2* tetha_j* exp_j *(dr[0]/drr2)*exp_i ) / (2 * m_sigmaAngularScale);
        }

        T force_x = -(dVLJ_dR * angularTerm + LJpot * dAngular_dx);
        T force_y = -(dVLJ_dR * angularTerm + LJpot * dAngular_dy);
        T torque = -(LJpot * dAngular_dPhi);
        // T torque = 0;

        return {{force_x, force_y, torque}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngularScale;
    T m_sigma6, m_sigma12;
};
