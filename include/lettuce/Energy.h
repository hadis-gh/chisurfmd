#pragma once

#include <vector>
#include <cmath>
#include <iomanip>

#include "lettuce/Vec.h"
#include "lettuce/constants.h"

template<typename TParticle, typename T = typename TParticle::value_type>
T calInternalTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T internalKE = calInternalKineticEnergy(particles, allSpecies);
    return (2 * internalKE) / (constants::boltzmann * particles.size() * 3.0);
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
    const T momentOfInertia = allSpecies[p.species].momentOfInertia;
    auto vel = getGeneralizedVelocities(p);

    for (size_t i = 0; i < vel.size(); ++i) {
        if (i < 2) {
            kineticEnergy += 0.5 * mass * vel[i] * vel[i];
        } else {
            kineticEnergy += 0.5 * momentOfInertia * vel[i] * vel[i];
        }
    }

    return kineticEnergy;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calParticleRelativeKineticEnergy(const TParticle& p, const std::vector<Species<T>>& allSpecies, const Vec<T>& relativeVelocity) {
    return 0.5 * allSpecies[p.species].mass * relativeVelocity.abs2();
}
