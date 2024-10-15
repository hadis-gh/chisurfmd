#pragma once

#include <vector>
#include <cmath>
#include <string>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/constants.h"
#include "lettuce/LennardJones.h"

template<typename T, template<typename TT> typename TParticle>
class AndersenThermostat {
public:
    AndersenThermostat(T dt, T collisionFrequency, T desiredTemperature, std::mt19937 gen)
        : dt(dt), collisionFrequency(collisionFrequency), desiredTemperature(desiredTemperature),
          dist(0.0, 1.0), maxwellDist(0.0, std::sqrt(constants::boltzmann * desiredTemperature)) {}

    void operator () (std::vector<TParticle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
        for (auto& p : particles) {
            if (dist(gen) < collisionFrequency * dt) {
                T mass = allSpecies[p.species].mass;
                p.v = maxwellDist(gen) / std::sqrt(mass);
            }
        }
    }

private:
    T dt;
    T collisionFrequency;
    T desiredTemperature;
    std::mt19937 gen;
    std::uniform_real_distribution<T> dist;
    std::normal_distribution<T> maxwellDist;
};

template<typename T, template<typename TT> typename TParticle>
class BerendsenThermostat {
public:
    BerendsenThermostat(T dt, T desiredTemperature, T relaxationTime)
        : dt(dt), desiredTemperature(desiredTemperature), relaxationTime(relaxationTime) {}

    void operator () (std::vector<TParticle<T>>& particles, const std::vector<Species<T>>& allSpecies) const {
        const T currentTemperature = calInternalTemperature(particles, allSpecies); 
        const auto lambda = std::sqrt(1 + dt / relaxationTime * (desiredTemperature / currentTemperature - 1));
        rescaleVelocity(particles, allSpecies, lambda);
    }

private:
    T dt;
    T relaxationTime;
    T desiredTemperature;
};

template<typename T, template<typename TT> typename TParticle>
class VelocityScalingThermostat {
public:
    VelocityScalingThermostat(T desiredTemperature)
        : desiredTemperature(desiredTemperature) {}

    void operator () (std::vector<TParticle<T>>& particles, const std::vector<Species<T>>& allSpecies) const {
        const T currentTemperature = calInternalTemperature(particles, allSpecies); 
        const auto lambda = std::sqrt(desiredTemperature / currentTemperature);
        rescaleVelocity(particles, allSpecies, lambda);
    }

private:
    T desiredTemperature;
};

template<typename T, typename TParticle>
T calInternalTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T internalKE = calInternalKineticEnergy(particles, allSpecies);
    return (2 * internalKE) / (constants::boltzmann * particles.size() * 3.0);
}

template<typename T, typename TParticle>
T calRawTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T rawKE = calRawKineticEnergy(particles, allSpecies);
    return (2 * rawKE) / (constants::boltzmann * particles.size() * 3.0);
}

template<typename T, typename TParticle>
T calInternalKineticEnergy(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comVel = calCOMVelocity(particles, allSpecies);
    T totalKE = 0.0;
    for (const auto& p : particles) {
        Vec<T> relativeVelocity = p.v - comVel;
        totalKE += calParticleRelativeKineticEnergy(p, allSpecies, relativeVelocity);
    }
    return totalKE;
}

template<typename T, typename TParticle>
T calRawKineticEnergy(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T totalKE = 0.0;
    for (const auto& p : particles) {
        totalKE += calParticleKineticEnergy(p, allSpecies);
    }
    return totalKE;
}

template<typename T, typename TParticle>
T calParticleKineticEnergy(const TParticle& p, const std::vector<Species<T>>& allSpecies) {
    return 0.5 * allSpecies[p.species].mass * p.v.abs2();
}

template<typename T, typename TParticle>
T calParticleRelativeKineticEnergy(const TParticle& p, const std::vector<Species<T>>& allSpecies, const Vec<T>& relativeVelocity) {
    return 0.5 * allSpecies[p.species].mass * relativeVelocity.abs2();
}

template<typename T, typename TParticle>
void rescaleVelocity(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T lambda) {
    for (auto& p : particles) {
        p.v *= lambda;
    }
}

template<typename T, typename TParticle>
Vec<T> calCOMVelocity(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> totalMomentum;
    T totalMass = 0.0;
    for (const auto& p : particles) {
        T mass = allSpecies[p.species].mass;
        totalMomentum += mass * p.v;
        totalMass += mass;
    }
    return totalMomentum / totalMass;
}

template<typename T, typename TParticle>
void removeCOMVelocity(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comVel = calCOMVelocity(particles, allSpecies);
    for (auto& p : particles) {
        p.v -= comVel;
    }
}