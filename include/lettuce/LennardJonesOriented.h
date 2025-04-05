#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/ParticleDot.h"
#include "lettuce/ParticleOriented.h"

template<typename TParticle, typename T = typename TParticle::value_type>
class LennardJonesOrientedForce {
public:
    using value_type = T;

    LennardJonesOrientedForce(T epsilon, T sigma, T cutoff, int phiOrder, T angularScale, T alpha)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6)
        , m_phiOrder(phiOrder), m_angularScale(angularScale)
        , m_alpha(alpha) {}

    Vec<T, 2> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {
        T const faceAng = std::atan2(dr[1], dr[0]);
        
        const T phi1new = p1.h * p1.d * (p1.phi - faceAng);
        const T phi2new = p2.h * p2.d * (p2.phi - faceAng);

        const T deltaPhi = phi2new - phi1new;

        if (r == 0 || r > m_cutoff) return {{0, 0}};
        
        const T min_distance = m_sigma * 0.5;
        const T effective_r = std::max(r, min_distance);

        const T r6 = effective_r * effective_r * effective_r * effective_r * effective_r * effective_r;
        const T r12 = r6 * r6;

        const T A = m_angularScale / r12;

        const T dA_dr = -A * 12 / effective_r;
        
        const T radialForce 
            = -4.0 * m_epsilon * (-12.0 * (m_sigma12 / r12) / effective_r + 6.0 * (m_sigma6 / r6) / effective_r) 
            - dA_dr * (1 + std::cos(m_phiOrder * deltaPhi + m_alpha));

        const T angularForce = m_phiOrder * A * std::sin(m_phiOrder * deltaPhi + m_alpha);
                
        return {{radialForce, angularForce}};
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
    int m_phiOrder;
    T m_angularScale;
    T m_alpha;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class LennardJonesOrientedPotential {
public:
    LennardJonesOrientedPotential(T epsilon, T sigma, T cutoff, int phiOrder, T angularScale, T alpha)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6)
        , m_phiOrder(phiOrder), m_angularScale(angularScale)
        , m_alpha(alpha) {}

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {
        T const faceAng = std::atan2(dr[1], dr[0]);

        const T phi1new = p1.h * p1.d * (p1.phi - faceAng);
        const T phi2new = p2.h * p2.d * (p2.phi - faceAng);

        const T deltaPhi = phi2new - phi1new;

        if (r == 0 || r > m_cutoff) return 0;

        const T r6 = r * r * r * r * r * r;
        const T r12 = r6 * r6;

        const T A = m_angularScale / r12;

        return 4.0 * m_epsilon * (m_sigma12 / r12 - m_sigma6 / r6) 
               + A * (1 + std::cos(m_phiOrder * deltaPhi + m_alpha));
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
    int m_phiOrder;
    T m_angularScale;
    T m_alpha;
};