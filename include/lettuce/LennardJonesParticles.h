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
#include "FieldCoupled.h"
#include "ChiralLJGeometric.h"
#include "TabularDFT.h"

namespace po = boost::program_options;

// ================================== Lennard-Jones V(r_ij) ,simple ==================================

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

// ================================== modified Lennard-Jones V(r, pi, pj) = V_LJ + cos(ΔΦ) ==================================

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

// ================================== Electric Field Idea V(r, pi, pj) = v_i*Q_i + v_j*Q_j ==================================

template<typename Particle, typename SFINAE = void>
struct FieldCoupled;

template<typename T>
struct FieldCoupled<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        LennardJones<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
            ("ChiralLJGeometricStrength", po::value<T>()->default_value(0.5), "Strength of chiral interaction");
    }

    static auto force(const po::variables_map &vm) {
        return FieldCoupledForce<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(), 
            vm["ChiralLJGeometricStrength"].as<T>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return FieldCoupledPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(), 
            vm["ChiralLJGeometricStrength"].as<T>()
        );
    }

    using ForceType = FieldCoupledForce<ParticleOriented<T>>;
};

// ================================== Lennard-Jones but replace r_ij with distance of amino acids V_LJ(R_ij) ... R_ij(r_ij, phi_i, phi_j) ==================================

template<typename Particle, typename SFINAE = void>
struct ChiralLJGeometric;

template<typename T>
struct ChiralLJGeometric<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        LennardJones<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
            ("ChiralLJGeometricStrength", po::value<T>()->default_value(0.5), "Strength of chiral interaction");
    }

    static auto force(const po::variables_map &vm) {
        return ChiralLJGeometricForce<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["exclusionRadius"].as<T>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return ChiralLJGeometricPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["exclusionRadius"].as<T>()
        );
    }

    using ForceType = ChiralLJGeometricForce<ParticleOriented<T>>;
};

// ================================== Extract Potential and Force from file ->> DFTB-based ==================================

template<typename Particle, typename SFINAE = void>
struct TabularDFT;

template<typename T>
struct TabularDFT<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        LennardJones<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
        ("tabularPotentialFile", po::value<std::string>() ,"file path of potential data")
        ("tabularForceFile", po::value<std::string>() ,"file path of force data")    
        ("bpIO", po::value<adios2::IO>() ,"adios IO")    
        ("parametersTabular", po::value<std::vector<T>>() ,"parameters for tabular potential")    
        ;
    }

    static auto force(const po::variables_map &vm) {
        return TabularDFTForce<ParticleOriented<T>>(
            vm["particlesType"].as<std::string>(),
            vm["tabularForceFile"].as<std::string>(),
            vm["bpIO"].as<adios2::IO>(),
            vm["parametersTabular"].as<std::vector<T>>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return TabularDFTPotential<ParticleOriented<T>>(
            vm["particlesType"].as<std::string>(),
            vm["tabularForceFile"].as<std::string>(),
            vm["bpIO"].as<adios2::IO>(),
            vm["parametersTabular"].as<std::vector<T>>()
        );
    }

    using ForceType = TabularDFT<ParticleOriented<T>>;
};