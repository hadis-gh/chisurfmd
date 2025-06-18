#pragma once

#include <vector>
#include <cmath>
#include "lettuce/core/Vec.h"
#include "lettuce/core/Circle.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

template<typename TParticle, typename T = typename TParticle::value_type>
class IsotropicLJForce {
public:
    using value_type = T;

    IsotropicLJForce(T epsilon, T sigma, T cutoff)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6) {}
    T operator()(const TParticle &p1, const TParticle &p2, const Vec<T, 2>& dr, const T r) const {
        if (r > m_cutoff) return 0;
        
        const T min_distance = m_sigma * 0.5;
        const T effective_r = std::max(r, min_distance);

        const T r6 = effective_r * effective_r * effective_r * effective_r * effective_r * effective_r;
        const T r12 = r6 * r6;

        return 48.0 * m_epsilon * (m_sigma12 / (r12 * effective_r) - 0.5 * m_sigma6 / (r6 * effective_r));
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class IsotropicLJPotential {
public:
    IsotropicLJPotential(T epsilon, T sigma, T cutoff)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6) {}

    T operator()(const TParticle &p1, const TParticle &p2, const Vec<T, 2>& dr, const T r) const {
        if (r == 0 || r > m_cutoff) return 0;

        const T r6 = r * r * r * r * r * r;
        const T r12 = r6 * r6;

        return 4.0 * m_epsilon * (m_sigma12 / r12 - m_sigma6 / r6);
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
};