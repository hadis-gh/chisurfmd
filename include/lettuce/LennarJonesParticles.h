#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/Particle.h"
#include "lettuce/ParticleOriented.h"
#include "LennardJones.h"
#include "LennardJonesOriented.h"

namespace po = boost::program_options;

template<typename Particle, typename SFINAE=void>
struct LennardJones;

template<typename T>
struct LennardJones<ParticleDot<T>>
{
    static void initProgramOptions(po::options_description &desc) 
    {
        desc.add_options()
            ("LJepsilon", po::value<T>()->default_value(1), "kkkk")
            ("LJsigma", po::value<T>()->default_value(1), "Phi for")
            ("LJcutoff", po::value<T>()->default_value(10.), "cutoff distance")
        ;
    }
    static auto force(const po::variables_map &vm)
    {
        return LennardJonesForce(vm["LJepsilon"].as<T>(), vm["LJsigma"].as<T>(), vm["LJcutoff"].as<T>());
    }
    
    static auto  potential(const po::variables_map &vm)
    {
        return LennardJonesPotential(vm["LJepsilon"].as<T>(), vm["LJsigma"].as<T>(), vm["LJcutoff"].as<T>());
    }

    using ForceType = LennardJonesForce<T>;
    // using ForceType = std::decay_t<decltype(force(std::declval<po::variables_map>()))>;
};

template<typename T>
struct LennardJones<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc)
    {
        LennardJones<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
            ("LJPhi", po::value<T>()->default_value(1), "Phi for")
        ;
    } 
    static auto force(const po::variables_map &vm)
    {
        return LennardJonesOrientedForce(vm["LJepsilon"].as<T>(), vm["LJsigma"].as<T>(), vm["LJcutoff"].as<T>(), vm["LJPhi"].as<T>());
    }
    static auto potential(const po::variables_map &vm)
    {
        return LennardJonesOrientedPotential(vm["LJepsilon"].as<T>(), vm["LJsigma"].as<T>(), vm["LJcutoff"].as<T>(), vm["LJPhi"].as<T>());
    }

    using ForceType = LennardJonesOrientedForce<T>;
};