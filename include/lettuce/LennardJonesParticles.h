#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/ParticleDot.h"
#include "lettuce/ParticleOriented.h"
#include "LennardJones.h"
#include "LennardJonesOriented.h"
#include "LennardJonesChiral.h"

namespace po = boost::program_options;

template<typename Particle, typename SFINAE = void>
struct LennardJones;

template<typename T>
struct LennardJones<ParticleDot<T>> 
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
            return LennardJonesForce<ParticleDot<T>>(
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
            return LennardJonesPotential<ParticleDot<T>>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>()
            );
        } catch (const boost::bad_any_cast& e) {
            std::cerr << "Error initializing LennardJonesPotential: " << e.what() << std::endl;
            throw;
        }
    }

    using ForceType = LennardJonesForce<ParticleDot<T>>;
};

template<typename T>
struct LennardJones<ParticleOriented<T>> 
{
    static void initProgramOptions(po::options_description &desc) {
        LennardJones<ParticleDot<T>>::initProgramOptions(desc); // Reuse base options
        desc.add_options()
            ("LJPhiOrder",      po::value<unsigned int>()   ->default_value(2),     "rotational order for orientation-dependent interactions")
            ("LJangularScale",  po::value<T>()              ->default_value(1.0),   "scaling factor for the orientation-dependent interaction")
            ("LJalpha",         po::value<T>()              ->default_value(M_PI),  "phase shift factor for the orientation-dependent interaction")
        ;
    }

    static auto force(const po::variables_map &vm) {
        try {
            return LennardJonesOrientedForce<ParticleOriented<T>>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>(), 
                vm["LJPhiOrder"].as<unsigned int>(),
                vm["LJangularScale"].as<T>(),
                vm["LJalpha"].as<T>()
            );
            std::cout << "scale of orientation: " << vm["LJangularScale"].as<T>() << std::endl;
            std::cout << "angular orientation order: " << vm["LJPhiOrder"].as<unsigned int>() << "\n\n";

        } catch (const boost::bad_any_cast& e) {
            std::cerr << "Error initializing LennardJonesOrientedForce: " << e.what() << std::endl;
            throw;
        }
    }

    static auto potential(const po::variables_map &vm) {
        try {
            return LennardJonesOrientedPotential<ParticleOriented<T>>(
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

    using ForceType = LennardJonesOrientedForce<ParticleOriented<T>>;
};

template<typename Particle, typename SFINAE = void>
struct LennardJonesChiral;

template<typename T>
struct LennardJonesChiral<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        LennardJones<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
            ("LJchiralStrength", po::value<T>()->default_value(0.5), "Strength of chiral interaction");
    }

    static auto force(const po::variables_map &vm) {
        return LennardJonesChiralForce<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(), 
            vm["LJchiralStrength"].as<T>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return LennardJonesChiralPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(), 
            vm["LJchiralStrength"].as<T>()
        );
    }

    using ForceType = LennardJonesChiralForce<ParticleOriented<T>>;
};