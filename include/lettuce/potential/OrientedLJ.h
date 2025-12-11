#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <boost/program_options.hpp>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"
#include "IsotropicLJ.h"

namespace po = boost::program_options;

template<typename TParticle, typename T = typename TParticle::value_type>
class OrientedLJPotential {
public:
    OrientedLJPotential(T epsilon, T sigma, T cutoff, int phiOrder, T angularScale, T alpha)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6)
        , m_phiOrder(phiOrder), m_angularScale(angularScale)
        , m_alpha(alpha) {}

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {
        
        const T deltaPhi = p2.phi - p1.phi;

        if (r == 0 || r > m_cutoff) return 0;

        const T min_distance = m_sigma * 0.5;
        const T effective_r = std::max(r, min_distance);

        const T r6 = effective_r * effective_r * effective_r * effective_r * effective_r * effective_r;
        const T r12 = r6 * r6;

        const T A = m_angularScale / r12;

        return 4.0 * m_epsilon * (m_sigma12 / r12 - m_sigma6 / r6) 
               + A * (1 + std::cos(m_phiOrder * deltaPhi + m_alpha));
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigma6, m_sigma12;
    int m_phiOrder;
    T m_angularScale, m_alpha;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class OrientedLJForce {
public:
    using value_type = T;

    OrientedLJForce(T epsilon, T sigma, T cutoff, int phiOrder, T angularScale, T alpha)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff)
        , m_sigma6(sigma * sigma * sigma * sigma * sigma * sigma)
        , m_sigma12(m_sigma6 * m_sigma6)
        , m_phiOrder(phiOrder), m_angularScale(angularScale)
        , m_alpha(alpha) {}

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {

        const T deltaPhi = p2.phi - p1.phi;

        if (r == 0 || r > m_cutoff) return {{0, 0, 0}};
        
        const T min_distance = m_sigma * 0.5;
        const T effective_r = std::max(r, min_distance);

        const T r6 = effective_r * effective_r * effective_r * effective_r * effective_r * effective_r;
        const T r12 = r6 * r6;

        const T A = m_angularScale / r12;

        const T dA_dr = -A * 12 / (effective_r);
        
        const T radialForce 
            = -4.0 * m_epsilon * (-12.0 * (m_sigma12 / r12) / effective_r + 6.0 * (m_sigma6 / r6) / effective_r) 
            - dA_dr * (1 + std::cos(m_phiOrder * deltaPhi + m_alpha));

        const T angularForce = m_phiOrder * A * std::sin(m_phiOrder * deltaPhi + m_alpha);
                
        return {{radialForce * dr[0] / r, radialForce * dr[1] / r, angularForce}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigma6, m_sigma12;
    int m_phiOrder;
    T m_angularScale, m_alpha;
};

// ================================== Factory Struct for OrientedLJ ==================================

template<typename Particle, typename SFINAE = void>
struct OrientedLJ;

template<typename T>
struct OrientedLJ<ParticleOriented<T>> 
{
    static void initProgramOptions(po::options_description &desc) {
        IsotropicLJ<ParticleDot<T>>::initProgramOptions(desc); // Reuse base options
        desc.add_options()
            ("LJPhiOrder",      po::value<unsigned int>()   ->default_value(2),     "rotational order for orientation-dependent interactions")
            ("LJangularScale",  po::value<T>()              ->default_value(1.0),   "scaling factor for the orientation-dependent interaction")
            ("LJalpha",         po::value<T>()              ->default_value(M_PI),  "phase shift factor for the orientation-dependent interaction")
        ;
    }

    static auto force(const po::variables_map &vm) {
        try {
            return OrientedLJForce<ParticleOriented<T>>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>(), 
                vm["LJPhiOrder"].as<unsigned int>(),
                vm["LJangularScale"].as<T>(),
                vm["LJalpha"].as<T>()
            );
        } catch (const boost::bad_any_cast& e) {
            std::cerr << "Error initializing LennardJonesOrientedForce: " << e.what() << std::endl;
            throw;
        }
    }

    static auto potential(const po::variables_map &vm) {
        try {
            return OrientedLJPotential<ParticleOriented<T>>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>(), 
                vm["LJPhiOrder"].as<unsigned int>(),
                vm["LJangularScale"].as<T>(),
                vm["LJalpha"].as<T>()
            );
        } catch (const boost::bad_any_cast& e) {
            std::cerr << "Error initializing LennardJonesOrientedPotential: " << e.what() << std::endl;
            throw;
        }
    }

    using ForceType = OrientedLJForce<ParticleOriented<T>>;
};