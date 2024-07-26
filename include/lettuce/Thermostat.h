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
    AndersenThermostat(T collisionFrequency, T desiredTemperature, std::mt19937 gen)
        : collisionFrequency(collisionFrequency), desiredTemperature(desiredTemperature),
          dist(0.0, 1.0), maxwellDist(0.0, std::sqrt(constants::boltzmann * desiredTemperature)) {}

    void operator () (std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, T dt) {
        for (auto& p : particles) {
            if (dist(gen) < collisionFrequency * dt) {
                T mass = allSpecies[p.species].mass;
                p.v *= maxwellDist(gen) / std::sqrt(mass);
            }
        }
    }

private:
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

    T operator () (const T &currentTemperature) const{
        return std::sqrt(1 + dt/ relaxationTime * (desiredTemperature/ currentTemperature - 1));
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

    T operator () (const T &currentTemperature) const{
        return std::sqrt(desiredTemperature/ currentTemperature);
    }

private:
    T desiredTemperature;
};

template<typename T>
T systemTemperature(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies){
    T totalKineticEnergy = systemKineticEnergy(particles, allSpecies);
    return (2 * totalKineticEnergy / (constants::boltzmann * particles.size() * 3.));
}

template<typename T>
T kineticEnergy(const Particle<T> &particle, const std::vector<Species<T>>& allSpecies){
    return 0.5 * allSpecies[particle.species].mass * particle.v.abs2();
}

template<typename T>
T systemKineticEnergy(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies){
    T totalKineticEnergy = 0.0;
    for (auto &p: particles){
        totalKineticEnergy += kineticEnergy(p, allSpecies);
    }
    return totalKineticEnergy;
}

template<typename T, typename ThermostatMethod>
void applyThermostat(std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const ThermostatMethod &thermostatMethod){
    T currentTemperature = systemTemperature(particles, allSpecies);
    if (currentTemperature <= 0) {
        std::cerr << "Error: Current temperature is non-positive: " << currentTemperature << std::endl;
        return;
    }    
    auto lambda = thermostatMethod(currentTemperature);
    for (auto &p : particles){
        p.v *= lambda;
    }
}

