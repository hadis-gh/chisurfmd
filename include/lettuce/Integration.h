#pragma once

#include <vector>
#include <cmath>
#include <string>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/LennardJones.h"
#include "lettuce/Thermostat.h"


template<typename T, typename Force>
void EulerStep(std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const T &dt, const Force &force){
    const auto accelerations = calAllAccelerations(particles, allSpecies, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        p.v += accelerations[i]*dt;
        p.r += p.v * dt;
    }    
}

template<typename T, typename Force>
void VelocityVerletStep(std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const T &dt, const Force &force){
    const auto old_accelerations = calAllAccelerations(particles, allSpecies, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        p.r += p.v * dt + old_accelerations[i]/2 * dt * dt;
    }
    const auto accelerations = calAllAccelerations(particles, allSpecies, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].v += (old_accelerations[i] + accelerations[i])/2 * dt;
    }
}    

template<typename T, typename Integrator, typename Force>
void integrate(std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, T dt, T Time, Force&& force, Integrator &&integrator) {
    int numSteps = static_cast<int>(Time / dt);

    for (size_t i = 0; i < numSteps; ++i){
        integrator(particles, allSpecies, dt, std::forward<Force>(force));
        thermostat(particles, allSpecies);
    }
}


template<typename T>
void writePositionToFile(const std::vector<Particle<T>> &particles, std::ostream &file) {
    for (const auto &p : particles) {
        file << p.r[0] << " " << p.r[1] << " ";
    }
    file << "\n";
}

template<typename T>
void writeKineticEToFile(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, std::ostream &file) {
    for (const auto &p : particles) {
        file << kineticEnergy(p, allSpecies) << " ";
    }
    file << "\n";
}

template<typename T, typename Potential>
void writePotentialEToFile(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, Potential &&potential, std::ostream &file){
    std::vector<Vec<T>> potentiaEnergy = calAllAccelerations(particles,allSpecies, std::forward<Potential> (potential));
    for (const auto &u : potentiaEnergy) {
        file << u.abs2() << " ";
    }
    file << "\n";
}