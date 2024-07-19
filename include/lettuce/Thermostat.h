#pragma once

#include <vector>
#include <cmath>
#include <string>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/constants.h"
#include "lettuce/LennardJones.h"



template<typename T>
class BerendensenThermostat {
public:
    BerendensenThermostat(T dt, T desiredTemperature, T relaxationTime)
        : dt(dt), desiredTemperature(desiredTemperature), relaxationTime(relaxationTime) {}

    T operator () (const T &currentTemperature) const{
        return sqrt(1 + dt/ relaxationTime * (desiredTemperature/ currentTemperature - 1));
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
        return sqrt(desiredTemperature/ currentTemperature);
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
    return 0.5 * allSpecies[particle.species].mass * particle.v.abs2() * particle.v.abs2();
}

template<typename T>
T systemKineticEnergy(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies){
    T totalKineticEnergy;
    for (auto &p: particles){
        totalKineticEnergy += kineticEnergy(p, allSpecies);
    }
    return totalKineticEnergy;
}

template<typename T, typename ThermostatMethod>
void applyThermostat(std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const ThermostatMethod &thermostatMethod){
    T currentTemperature = systemTemperature(particles, allSpecies);
    auto lamda = thermostatMethod(currentTemperature);
    for (auto &p : particles){
        p.v *= lamda;
    }
}

