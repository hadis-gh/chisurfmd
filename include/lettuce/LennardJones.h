#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/ParticleDot.h"

template<typename T>
class LennardJonesForce {
public:
    using value_type = T;

    LennardJonesForce(T epsilon, T sigma, T cutoff)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6) {}

    T operator()(const T r) const {
        if (r == 0 || r > m_cutoff) return 0;

        const T r6 = r * r * r * r * r * r;
        const T r12 = r6 * r6;

        return 48.0 * m_epsilon * (m_sigma12 / (r12 * r) - 0.5 * m_sigma6 / (r6 * r));
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
};

template<typename T>
class LennardJonesPotential {
public:
    LennardJonesPotential(T epsilon, T sigma, T cutoff)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6) {}

    T operator()(const T r) const {
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