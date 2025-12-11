#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <boost/program_options.hpp>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"
#include "IsotropicLJ.h"

// Include adios2 if needed
#include <adios2.h>

namespace po = boost::program_options;

template<typename TParticle, typename T = typename TParticle::value_type>
class TabularDFTForce {
public:
    using value_type = T;

    TabularDFTForce(const std::string& particlesType, const std::string& adiosInput, adios2::IO& bpIO, const std::vector<T>& params)
        : m_particlesType(particlesType), m_adiosInput(adiosInput), m_bpIO(bpIO), m_params(params) {}

    Vec<T, 2> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {
        T const faceAng = std::atan2(dr[1], dr[0]);
        
        const T phi1new = p1.h * p1.d * (p1.phi - faceAng); 
        const T phi2new = p2.h * p2.d * (p2.phi - faceAng); 

        const T deltaPhi = phi2new - phi1new;
        
        const T radialForce = 0;

        const T angularForce = 0;
                
        return {{radialForce, angularForce}};
    }

private:
    std::string m_particlesType;
    std::string m_adiosInput;
    adios2::IO m_bpIO;
    std::vector<T> m_params;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class TabularDFTPotential {
public:
    TabularDFTPotential(const std::string& particlesType, const std::string& adiosInput, adios2::IO& bpIO, const std::vector<T>& params)
    : m_particlesType(particlesType), m_adiosInput(adiosInput), m_bpIO(bpIO), m_params(params) {}

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {
        T const faceAng = std::atan2(dr[1], dr[0]);

        const T phi1new = p1.h * p1.d * (p1.phi - faceAng); 
        const T phi2new = p2.h * p2.d * (p2.phi - faceAng); 

        return 0.0;
    }

    private:
    std::string m_particlesType;
    std::string m_adiosInput;
    adios2::IO m_bpIO;
    std::vector<T> m_params;
};

// ================================== Factory Struct for TabularDFT ==================================

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

    using ForceType = TabularDFTForce<ParticleOriented<T>>;
};