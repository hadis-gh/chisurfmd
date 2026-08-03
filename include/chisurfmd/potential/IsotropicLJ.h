#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <boost/program_options.hpp>
#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/Circle.h"
#include "chisurfmd/core/ParticleDot.h"
#include "chisurfmd/core/ParticleOriented.h"

namespace po = boost::program_options;

template<typename T>
inline T lennardJones(const T& r, const T& sigma, const T& epsilon) {
    if (r <= 0) return std::numeric_limits<T>::infinity();
    const T sr = sigma / r;
    const T sr6 = sr * sr * sr * sr * sr * sr;
    const T sr12 = sr6 * sr6;
    return 4.0 * epsilon * (sr12 - sr6);
}

template<typename T>
inline T lennardJonesDerivative(const T& r, const T& sigma, const T& epsilon) {
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
    
    Vec<T, 2> operator()(const TParticle &p1, const TParticle &p2, const Vec<T, 2>& dr, const T r) const {
        if (r > m_cutoff) return {{0, 0}};

        T force_magnitude = -lennardJonesDerivative(r, m_sigma, m_epsilon);
        return {{force_magnitude * dr[0] / r, force_magnitude * dr[1] / r}};
    }

private:
    T m_epsilon;
    T m_sigma;
    T m_cutoff;
    T m_sigma6;
    T m_sigma12;
};

// ================================== Factory Struct for IsotropicLJ ==================================

template<typename Particle, typename SFINAE = void>
struct IsotropicLJ;

template<typename T>
struct IsotropicLJ<ParticleDot<T>> 
{
    static void initProgramOptions(po::options_description &desc) {
        desc.add_options()
            ("LJepsilon", po::value<T>()->default_value(1.),     "epsilon in Lennard-Jones force and potential")
            ("LJsigma",   po::value<T>()->default_value(1.),     "sigma in Lennard-Jones force and potential")
            ("LJcutoff",  po::value<T>()->default_value(10.),    "cutoff distance for Lennard-Jones interactions")
        ;
    }

    static auto force(const po::variables_map &vm) {
        try {
            return IsotropicLJForce<ParticleDot<T>>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>()
            );
        } catch (const boost::bad_any_cast& e) {
            std::cerr << "Error initializing LennardJonesForce: " << e.what() << std::endl;
            throw;
        }
    }

    static auto potential(const po::variables_map &vm) {
        try {
            return IsotropicLJPotential<ParticleDot<T>>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>()
            );
        } catch (const boost::bad_any_cast& e) {
            std::cerr << "Error initializing LennardJonesPotential: " << e.what() << std::endl;
            throw;
        }
    }

    using ForceType = IsotropicLJForce<ParticleDot<T>>;
};