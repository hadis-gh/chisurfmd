#pragma once

#include <vector>
#include <cmath>
#include <string>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/LennardJones.h"


template<typename T>
T systemTemperature(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies){
    T totalKineticEnergy = systemKineticEnergy(particles, allSpecies);
    T boltzmanC = 1;
    return (2 * totalKineticEnergy / (boltzmanC * particles.size() * 3.));
}

template<typename T>
struct systemTemprature {
    systemTemprature(T boltzmanConstant) 
        : boltzmanConstant(boltzmanConstant) {}

    T operator()(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies){
        T totalKineticEnergy = systemKineticEnergy(particles, allSpecies);
        return (2 * totalKineticEnergy / (boltzmanConstant * particles.size() * 3.));
    }
private:
    T boltzmanConstant;
};

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

template<typename T>
T berendensonT (){

    currentTemperature = systemKineticEnergy();
    return (1 + dt/ relaxationT * (desireTemperature/ currentTemperature - 1));
}

template<typename T, typename Temperature, typename ThermostatMethod>
void thermostat(std::vector<Particle<T>> &particles, const ThermostatMethod &thermostatMethod){
    auto lamda = thermostatMethod;
    for (auto &p : particles){
        p.v *= lamda;
    }
}

