#pragma once

#include <vector>
#include <cmath>
#include "lettuce/core/Vec.h"
#include "lettuce/core/Circle.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

template<typename T>
constexpr T lennardJones(const T& r, const T& sigma, const T& epsilon) {
    if (r <= 0) return std::numeric_limits<T>::infinity();
    const T sr = sigma / r;
    const T sr6 = sr * sr * sr * sr * sr * sr;
    const T sr12 = sr6 * sr6;
    return 4.0 * epsilon * (sr12 - sr6);
}

template<typename T>
constexpr T lennardJonesDerivative(const T& r, const T& sigma, const T& epsilon) {
    if (r <= 0) return std::numeric_limits<T>::infinity();
    const T sr = sigma / r;
    const T sr6 = sr * sr * sr * sr * sr * sr;
    const T sr12 = sr6 * sr6;
    return 4.0 * epsilon * (6.0 * sr6 / r - 12.0 * sr12 / r);
}

template<typename TParticle, typename T = typename TParticle::value_type>
class IsotropicLJPotential {
public:
    IsotropicLJPotential(T epsilon, T sigma, T cutoff)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6) {}

    T operator()(const TParticle &p1, const TParticle &p2, const Vec<T, 2>& dr, const T r) const {
        if (r == 0 || r > m_cutoff) return 0;

        return lennardJones(r, m_sigma, m_epsilon);
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
};

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

        return -lennardJonesDerivative(r, m_sigma, m_epsilon) * dr/r;
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
};