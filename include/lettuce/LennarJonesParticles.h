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
            ("LJcutoff",  po::value<T>()->default_value(10.),   "cutoff distance for Lennard-Jones interactions")
        ;
    }

    static auto force(const po::variables_map &vm) {
        try {
            return LennardJonesForce<T>(
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
            return LennardJonesPotential<T>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>()
            );
        } catch (const boost::bad_any_cast& e) {
            std::cerr << "Error initializing LennardJonesPotential: " << e.what() << std::endl;
            throw;
        }
    }

    using ForceType = LennardJonesForce<T>;
};

template<typename T>
struct LennardJones<ParticleOriented<T>> 
{
    static void initProgramOptions(po::options_description &desc) {
        LennardJones<ParticleDot<T>>::initProgramOptions(desc); // Reuse base options
        desc.add_options()
            ("LJPhiOrder", po::value<int>()->default_value(4), "rotational order for orientation-dependent interactions")
            ("LJangularScale", po::value<T>()->default_value(5.0), "scaling factor for the orientation-dependent interaction")
        ;
    }

    static auto force(const po::variables_map &vm) {
        try {
            return LennardJonesOrientedForce<T>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>(), 
                vm["LJPhiOrder"].as<int>(),
                vm["LJangularScale"].as<T>()
            );
            std::cout << "scale of orientation: " << vm["LJangularScale"].as<T>() << std::endl;
            std::cout << "angular orientation order: " << vm["LJPhiOrder"].as<int>() << "\n\n";

        } catch (const boost::bad_any_cast& e) {
            std::cerr << "Error initializing LennardJonesOrientedForce: " << e.what() << std::endl;
            throw;
        }
    }

    static auto potential(const po::variables_map &vm) {
        try {
            return LennardJonesOrientedPotential<T>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>(), 
                vm["LJPhiOrder"].as<int>(),
                vm["LJangularScale"].as<T>()
            );
        } catch (const boost::bad_any_cast& e) {
            std::cerr << "Error initializing LennardJonesOrientedPotential: " << e.what() << std::endl;
            throw;
        }
    }

    using ForceType = LennardJonesOrientedForce<T>;
};
