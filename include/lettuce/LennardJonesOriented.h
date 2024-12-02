#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/ParticleDot.h"

template<typename T>
class LennardJonesOrientedForce {
public:
    using value_type = T;

    LennardJonesOrientedForce(T epsilon, T sigma, T cutoff, int phiOrder)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6)
        , m_phiOrder(phiOrder) {}

    Vec<T, 2> operator()(const T r, const T deltaPhi) const {
        if (r == 0 || r > m_cutoff) return {{0, 0}};

        const T r6 = r * r * r * r * r * r;
        const T r12 = r6 * r6;

        const T A = 5 / r12;

        const T dA_dr = -A * 6 / r;

        const T radialForce 
            = -4.0 * m_epsilon * (-12.0 * (m_sigma12 / r12) / r + 6.0 * (m_sigma6 / r6) / r) 
            + dA_dr * std::cos(m_phiOrder * deltaPhi);

        const T angularForce = m_phiOrder * A * std::sin(m_phiOrder * deltaPhi);
        
        return {{radialForce, angularForce}};
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
    int m_phiOrder;
};

template<typename T>
class LennardJonesOrientedPotential {
public:
    LennardJonesOrientedPotential(T epsilon, T sigma, T cutoff, int phiOrder)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6)
        , m_phiOrder(phiOrder) {}

    T operator()(const T r, const T deltaPhi) const {
        if (r == 0 || r > m_cutoff) return 0;

        const T r6 = r * r * r * r * r * r;
        const T r12 = r6 * r6;

        constexpr T A_factor = 5;
        const T A = A_factor / r12;

        return 4.0 * m_epsilon * (m_sigma12 / r12 - m_sigma6 / r6) 
               + A * std::cos(m_phiOrder * deltaPhi);
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
    int m_phiOrder;
};