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
#include "TabularPotential.h"

namespace po = boost::program_options;

template<typename Particle, typename SFINAE = void>
struct PotentialFactory;

template<typename T>
struct PotentialFactory<ParticleDot<T>> 
{
    static void initProgramOptions(po::options_description &desc) {
        desc.add_options()
            ("potentialType",         po::value<std::string>()->default_value("LennardJones"),    "Type of potential (LennardJones, TabularPotential)")
            ("LJepsilon",             po::value<T>()->default_value(1.),                          "epsilon in Lennard-Jones force and potential")
            ("LJsigma",               po::value<T>()->default_value(1.),                          "sigma in Lennard-Jones force and potential")
            ("LJcutoff",              po::value<T>()->default_value(10.),                         "cutoff distance for Lennard-Jones interactions")
            ("TabularPairType",       po::value<std::string>()->default_value("RR_UU"),           "Pair type for tabular potential")
            ("TabularPFile",          po::value<std::string>()->default_value("potential.bp"),    "File containing tabular potential data")
            ("TabularPParams",        po::value<std::vector<T>>()->multitoken()->default_value(std::vector<T>{NAN, NAN, 0., NAN}), "Parameters for tabular potential")
        ;
    }

    static auto force(const po::variables_map &vm) {
        std::string potentialType = vm["potentialType"].as<std::string>();
        if (potentialType == "LennardJones") {
            return LennardJonesForce<T>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>()
            );
        } else if (potentialType == "TabularPotential") {
            adios2::ADIOS adios;
            adios2::IO bpIO = adios.DeclareIO("potential");
            return TabularPotentialForce<T>(
                vm["TabularPairType"].as<std::string>(), 
                vm["TabularPotentialFile"].as<std::string>(), 
                bpIO, 
                vm["TabularPotentialParams"].as<std::vector<T>>()
            );
        } else {
            throw std::runtime_error("Unknown potential type");
        }
    }

    static auto potential(const po::variables_map &vm) {
        std::string potentialType = vm["potentialType"].as<std::string>();
        if (potentialType == "LennardJones") {
            return LennardJonesPotential<T>(
                vm["LJepsilon"].as<T>(), 
                vm["LJsigma"].as<T>(), 
                vm["LJcutoff"].as<T>()
            );
        } else if (potentialType == "TabularPotential") {
            adios2::ADIOS adios;
            adios2::IO bpIO = adios.DeclareIO("potential");
            return TabularPotentialPotential<T>(
                vm["TabularPairType"].as<std::string>(), 
                vm["TabularPotentialFile"].as<std::string>(), 
                bpIO, 
                vm["TabularPotentialParams"].as<std::vector<T>>()
            );
        } else {
            throw std::runtime_error("Unknown potential type");
        }
    }

    using ForceType = LennardJonesOrientedForce<T>;
};