#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/ParticleDot.h"
#include "lettuce/ParticleOriented.h"
#include "LennardJones.h"
#include "LennardJonesOriented.h"

namespace po = boost::program_options;

template<typename Particle, typename SFINAE=void>
struct LennardJones;

template<typename T>
struct LennardJones<ParticleDot<T>> 
{
    static void initProgramOptions(po::options_description &desc) {
        desc.add_options()
            ("LJepsilon", po::value<T>()->default_value(1),     "epsilon in Lennard-Jones force and potential")
            ("LJsigma",   po::value<T>()->default_value(1),     "sigma in Lennard-Jones force and potential")
            ("LJcutoff",  po::value<T>()->default_value(10.),   "cutoff distance for Lennard-Jones interactions");
    }

    static auto force(const po::variables_map &vm) {
        return LennardJonesForce<T>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return LennardJonesPotential<T>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>()
        );
    }

    using ForceType = LennardJonesForce<T>;
};

template<typename T>
struct LennardJones<ParticleOriented<T>> 
{
    static void initProgramOptions(po::options_description &desc) {
        LennardJones<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
            ("LJPhiOrder", po::value<int>()->default_value(4), "rotational order for orientation-dependent interactions");
    }

    static auto force(const po::variables_map &vm) {
        return LennardJonesOrientedForce<T>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(), 
            vm["LJPhiOrder"].as<int>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return LennardJonesOrientedPotential<T>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(), 
            vm["LJPhiOrder"].as<int>()
        );
    }

    using ForceType = LennardJonesOrientedForce<T>;
};