#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/ParticleDot.h"
#include "lettuce/ParticleOriented.h"

template<typename TParticle, typename T = typename TParticle::value_type>
class FieldCoupledForce {
public:
    FieldCoupledForce(T epsilon, T sigma, T cutoff, T chiralStrength)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_chiralStrength(chiralStrength)
        , m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)) {}

    Vec<T, 2> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {
        if (r == 0 || r > m_cutoff) return {{0, 0}};

        T r6 = std::pow(r, 6);
        T r12 = r6 * r6;
        T LJforce = -4 * m_epsilon * (-12.0 * m_sigma12 / r12 + 6.0 * m_sigma6 / r6) / r;

        T gamma = std::atan2(dr[1], dr[0]);

        T psi_i = p1.phi - gamma;
        T psi_j = p2.phi - gamma;

        T Qi = std::sin(psi_i);
        T Qj = std::sin(psi_j);

        T Fi = m_chiralStrength / (r * r);
        T Fj = Fi; // symmetric ? asymmetric later

        T vi = Fi * Qi;
        T vj = Fj * Qj;

        T V = vi * Qj + vj * Qi;

        // Derivative of V fix later
        T chiralForce = -V / r;

        T totalForce = LJforce + chiralForce;

        return {{LJforce, chiralForce}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff;
    T m_sigma6, m_sigma12;
    T m_chiralStrength;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class FieldCoupledPotential {
public:
    FieldCoupledPotential(T epsilon, T sigma, T cutoff, T chiralStrength)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_chiralStrength(chiralStrength)
        , m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)) {}

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {
        if (r == 0 || r > m_cutoff) return 0;

        T r6 = std::pow(r, 6);
        T r12 = r6 * r6;
        T LJ = 4 * m_epsilon * (m_sigma12 / r12 - m_sigma6 / r6);

        T gamma = std::atan2(dr[1], dr[0]);
        T psi_i = p1.phi - gamma;
        T psi_j = p2.phi - gamma;

        T Qi = std::sin(psi_i);
        T Qj = std::sin(psi_j);
        T Fi = m_chiralStrength / (r * r);
        T vi = Fi * Qi;
        T vj = Fi * Qj;

        T V = vi * Qj + vj * Qi;

        return LJ + V;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_chiralStrength;
    T m_sigma6, m_sigma12;
};