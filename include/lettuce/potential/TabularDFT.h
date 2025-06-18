#pragma once

#include <vector>
#include <cmath>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

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