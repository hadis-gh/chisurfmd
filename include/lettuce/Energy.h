#pragma once

#include <vector>
#include <cmath>
#include <iomanip>

#include "lettuce/Vec.h"
#include "lettuce/constants.h"

// template<typename TParticle, typename T = typename TParticle::value_type>
// T calInternalTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
//     T internalKE = calInternalKineticEnergy(particles, allSpecies);
//     return (2 * internalKE) / (constants::boltzmann * particles.size() * 3.0);
// }

template<typename TParticle, typename T = typename TParticle::value_type>
T calInternalTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T translationalKE = 0.0;
    T rotationalKE = 0.0;
    size_t numParticles = particles.size();
    
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        const T MI = allSpecies[p.species].momentOfInertia;
        auto vel = getGeneralizedVelocities(p);

        for (size_t i = 0; i < 2; ++i) {
            translationalKE += 0.5 * mass * vel[i] * vel[i];
        }

        if (vel.size() > 2) {
            rotationalKE += 0.5 * MI * vel[2] * vel[2];
        }
    }

    // Each degree of freedom contributes 1/2 kB T to the energy: equipartition theorem
    T totalDOF = 2 * numParticles;  // 2 translational DOF per particle
    if (DegreesOfFreedom<TParticle>::degreesOfFreedom() > 2) {
        totalDOF += numParticles;    // Add 1 rotational DOF per particle if present
    }

    return 2.0 * (translationalKE + rotationalKE) / (constants::boltzmann * totalDOF);
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calRawTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T rawKE = calRawKineticEnergy(particles, allSpecies);
    return (2 * rawKE) / (constants::boltzmann * particles.size() * 3.0);
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calInternalKineticEnergy(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comVel = calCOMVelocity(particles, allSpecies);
    T totalKE = 0.0;
    for (const auto& p : particles) {
        Vec<T> relativeVelocity = p.v - comVel;
        totalKE += calParticleRelativeKineticEnergy(p, allSpecies, relativeVelocity);
    }
    return totalKE;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calRawKineticEnergy(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T totalKE = 0.0;
    for (const auto& p : particles) {
        totalKE += calParticleKineticEnergy(p, allSpecies);
    }
    return totalKE;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calParticleKineticEnergy(const TParticle& p, const std::vector<Species<T>>& allSpecies) {
    T kineticEnergy = 0.0;
    const T mass = allSpecies[p.species].mass;
    const T MI = allSpecies[p.species].momentOfInertia;
    auto vel = getGeneralizedVelocities(p);

    for (size_t i = 0; i < vel.size(); ++i) {
        if (i < 2) {
            kineticEnergy += 0.5 * mass * vel[i] * vel[i];
        } else {
            kineticEnergy += 0.5 * MI * vel[i] * vel[i];
        }
    }
    return kineticEnergy;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calParticleRelativeKineticEnergy(const TParticle& p, const std::vector<Species<T>>& allSpecies, const Vec<T>& relativeVelocity) {
    return 0.5 * allSpecies[p.species].mass * relativeVelocity.abs2();
}
