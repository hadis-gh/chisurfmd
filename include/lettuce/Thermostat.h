#pragma once

#include <vector>
#include <cmath>
#include <string>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/constants.h"
#include "lettuce/LennardJones.h"

template<typename T>
class AndersenThermostat {
public:
    AndersenThermostat(T dt, T collisionFrequency, T desiredTemperature, std::mt19937 gen)
        : dt(dt), collisionFrequency(collisionFrequency), desiredTemperature(desiredTemperature),
          dist(0.0, 1.0), maxwellDist(0.0, std::sqrt(constants::boltzmann * desiredTemperature)) {}

    void operator () (std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
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

template<typename T>
class BerendsenThermostat {
public:
    BerendsenThermostat(T dt, T desiredTemperature, T relaxationTime)
        : dt(dt), desiredTemperature(desiredTemperature), relaxationTime(relaxationTime) {}

    void operator () (std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) const {
        const T currentTemperature = calculateInternalTemperature(particles, allSpecies); 
        const auto lambda = std::sqrt(1 + dt / relaxationTime * (desiredTemperature / currentTemperature - 1));
        rescaleVelocity(particles, allSpecies, lambda);
    }

private:
    T dt;
    T relaxationTime;
    T desiredTemperature;
};

template<typename T>
class VelocityScalingThermostat {
public:
    VelocityScalingThermostat(T desiredTemperature)
        : desiredTemperature(desiredTemperature) {}

    void operator () (std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) const {
        const T currentTemperature = calculateInternalTemperature(particles, allSpecies); 
        const auto lambda = std::sqrt(desiredTemperature / currentTemperature);
        rescaleVelocity(particles, allSpecies, lambda);
    }

private:
    T desiredTemperature;
};

template<typename T>
T calculateInternalTemperature(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
    T internalKE = calculateInternalKineticEnergy(particles, allSpecies);
    return (2 * internalKE) / (constants::boltzmann * particles.size() * 3.0);
}

template<typename T>
T calculateRawTemperature(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
    T rawKE = calculateRawKineticEnergy(particles, allSpecies);
    return (2 * rawKE) / (constants::boltzmann * particles.size() * 3.0);
}

template<typename T>
T calculateInternalKineticEnergy(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comVel = calculateCOMVelocity(particles, allSpecies);
    T totalKE = 0.0;
    for (const auto& p : particles) {
        Vec<T> relativeVelocity = p.v - comVel;
        totalKE += calculateParticleRelativeKineticEnergy(p, allSpecies, relativeVelocity);
    }
    return totalKE;
}

template<typename T>
T calculateRawKineticEnergy(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
    T totalKE = 0.0;
    for (const auto& p : particles) {
        totalKE += calculateParticleKineticEnergy(p, allSpecies);
    }
    return totalKE;
}

template<typename T>
T calculateParticleKineticEnergy(const Particle<T>& p, const std::vector<Species<T>>& allSpecies) {
    return 0.5 * allSpecies[p.species].mass * p.v.abs2();
}

template<typename T>
T calculateParticleRelativeKineticEnergy(const Particle<T>& p, const std::vector<Species<T>>& allSpecies, const Vec<T>& relativeVelocity) {
    return 0.5 * allSpecies[p.species].mass * relativeVelocity.abs2();
}

template<typename T>
void rescaleVelocity(std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, const T lambda) {
    for (auto& p : particles) {
        p.v *= lambda;
    }
}

template<typename T>
Vec<T> calculateCOMVelocity(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> totalMomentum;
    T totalMass = 0.0;
    for (const auto& p : particles) {
        T mass = allSpecies[p.species].mass;
        totalMomentum += mass * p.v;
        totalMass += mass;
    }
    return totalMomentum / totalMass;
}

template<typename T>
void removeCOMVelocity(std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comVel = calculateCOMVelocity(particles, allSpecies);
    for (auto& p : particles) {
        p.v -= comVel;
    }
}
