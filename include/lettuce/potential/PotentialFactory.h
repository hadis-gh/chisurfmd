#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include "lettuce/core/Vec.h"
#include "lettuce/core/Circle.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"
#include "IsotropicLJ.h"
#include "OrientedLJ.h"
#include "FieldCoupled.h"
#include "GeometricLJ.h"
#include "TabularDFT.h"
#include "PatchyLJ.h"

namespace po = boost::program_options;

// ================================== Lennard-Jones V(r_ij) ,simple ==================================

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

// ================================== modified Lennard-Jones V(r, pi, pj) = V_LJ + cos(ΔΦ) ==================================

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
            std::cout << "scale of orientation: " << vm["LJangularScale"].as<T>() << std::endl;
            std::cout << "angular orientation order: " << vm["LJPhiOrder"].as<unsigned int>() << "\n\n";

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

// ================================== Electric Field Idea V(r, pi, pj) = v_i*Q_i + v_j*Q_j ==================================

template<typename Particle, typename SFINAE = void>
struct FieldCoupled;

template<typename T>
struct FieldCoupled<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        IsotropicLJ<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
            ("fieldCoupleFactor", po::value<T>()->default_value(0.5), "Strength of chiral interaction");
    }

    static auto force(const po::variables_map &vm) {
        return FieldCoupledForce<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(), 
            vm["fieldCoupleFactor"].as<T>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return FieldCoupledPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(), 
            vm["fieldCoupleFactor"].as<T>()
        );
    }

    using ForceType = FieldCoupledForce<ParticleOriented<T>>;
};

// ================================== Lennard-Jones but replace r_ij with distance of amino acids V_LJ(R_ij) ... R_ij(r_ij, phi_i, phi_j) ==================================

template<typename Particle, typename SFINAE = void>
struct GeometricLJ;

template<typename T>
struct GeometricLJ<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        IsotropicLJ<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
            ("sigmaPatchyLJ", po::value<T>()->default_value(0.5), "Strength of chiral interaction");
    }

    static auto force(const po::variables_map &vm) {
        return GeometricLJForce<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["exclusionRadius"].as<T>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return GeometricLJPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["exclusionRadius"].as<T>()
        );
    }

    using ForceType = GeometricLJForce<ParticleOriented<T>>;
};

// ================================== Patchy colloidal potential - multiply of Lennard-Jones with angular term ==================================

template<typename Particle, typename SFINAE = void>
struct PatchyLJ;

template<typename T>
struct PatchyLJ<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        IsotropicLJ<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
            ("sigmaAngularScale",   po::value<T>()->default_value(4.),              "asitropic strength of potential")
            ("patchNum",            po::value<unsigned int>()->default_value(3),    "asitropic strength of potential")
            ("useSumPatch",         po::value<bool>()->default_value(false),        "asitropic strength of potential")    
        ;
    }

    static auto force(const po::variables_map &vm) {
        return PatchyLJForce<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["sigmaAngularScale"].as<T>(),
            vm["patchNum"].as<unsigned int>(),
            vm["useSumPatch"].as<bool>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return PatchyLJPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["sigmaAngularScale"].as<T>(),
            vm["patchNum"].as<unsigned int>(),
            vm["useSumPatch"].as<bool>()
        );
    }

    using ForceType = PatchyLJForce<ParticleOriented<T>>;
};

// ================================== Extract Potential and Force from file ->> DFTB-based ==================================

template<typename Particle, typename SFINAE = void>
struct TabularDFT;

template<typename T>
struct TabularDFT<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        IsotropicLJ<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
        ("tabularPotentialFile",    po::value<std::string>(),      "file path of potential data")
        ("tabularForceFile",        po::value<std::string>(),      "file path of force data")    
        ("bpIO",                    po::value<adios2::IO>(),       "adios IO")    
        ("parametersTabular",       po::value<std::vector<T>>(),   "parameters for tabular potential")    
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