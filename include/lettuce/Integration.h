#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/LennardJones.h"
#include "lettuce/Thermostat.h"

template<typename T>
void implementPBC(Particle<T> &p, const T& boxPBC){
    if (p.r[0] > boxPBC) { p.r[0] -= boxPBC; }
    else if (p.r[0] < 0) { p.r[0] += boxPBC; }
    
    if (p.r[1] > boxPBC) { p.r[1] -= boxPBC; }
    else if (p.r[1] < 0) { p.r[1] += boxPBC; }
}

template<typename T, typename Force>
void EulerSymplecticStep(std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const T &dt, const T& boxPBC, const Force &force){
    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        p.v += accelerations[i]*dt;
        p.r += p.v * dt;
        implementPBC(p, boxPBC);
    }    
}

template<typename T, typename Force>
void EulerStep(std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const T &dt, const T& boxPBC, const Force &force){
    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        p.r += p.v * dt;
        p.v += accelerations[i]*dt;
        implementPBC(p, boxPBC);
    }    
}

template<typename T, typename Force>
void VelocityVerletStep(std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const T &dt, const T& boxPBC, const Force &force) {
    const auto old_accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    std::vector<T> torques(particles.size(), 0);
    
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        p.r += p.v * dt + old_accelerations[i]/2 * dt * dt;
        implementPBC(p, boxPBC);
    }

    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].v += (old_accelerations[i] + accelerations[i]) / 2 * dt;
    }
}    

template<typename T, typename Integrator, typename Force>
void integrate(std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, T dt, T Time, const T&boxPBC,
                Force&& force, Integrator &&integrator) {
    int numSteps = static_cast<int>(Time / dt);

    for (size_t i = 0; i < numSteps; ++i){
        integrator(particles, allSpecies, dt, boxPBC, std::forward<Force>(force));
    }
}

template<typename T>
T averageNeighbors(const std::vector<Particle<T>> &particles, const T &distance){
    T totalNeighbors = 0;
    for (auto p1: particles){
        int eachParticleNeighbors = 0;
        for (auto p2: particles){
            if (p1.r != p2.r && ((p1.r - p2.r).abs()) <= distance){
                ++ eachParticleNeighbors;
            }
        }
        totalNeighbors += eachParticleNeighbors;
    }
    return totalNeighbors/ particles.size();
}

template<typename T>
void writeInitialParticles(const std::vector<Particle<T>>& particles, const T radius) {
    std::ofstream initialParticles("FirsConfigPlot.dat");
    for (auto &p: particles){
        initialParticles << p.r[0] << " " << p.r[1] << " " << radius << std::endl;
    }
}

template<typename T>
void writePositionToFile(const std::vector<Particle<T>> &particles, std::ostream &file, const T &dt, const int &step) {
    file << dt * step << " ";

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

template<typename T>
void writeRelativeKineticEToFile(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, std::ostream &file) {
    Vec<T> comVelocity = computeCenterOfMassVelocity(particles, allSpecies);

    for (const auto &p : particles) {
        Vec<T> relativeVelocity = p.v - comVelocity;
        file << relativeKineticEnergy(p, allSpecies, relativeVelocity) << " ";
    }
    file << "\n";
}

template<typename T>
void writeTemperature(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, std::ostream &file) {
    file << systemTemperature(particles, allSpecies) << std::endl;
}

template<typename T>
void writeRealTemperature(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, std::ostream &file) {
    file << systemTemperatureOld(particles, allSpecies) << std::endl;
}

template<typename T, typename Potential>
void writePotentialEToFile(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Potential &&potential, std::ostream &file){
    std::vector<Vec<T>> potentiaEnergy = calAllAccelerations(particles,allSpecies, boxPBC, std::forward<Potential> (potential));
    for (const auto &u : potentiaEnergy) {
        file << u.abs() << " ";
    }
    file << "\n";
}

template<typename T>
void writeAverageNeighborToFile(const std::vector<Particle<T>> &particles, const  std::vector<T> &distances, std::ostream &file, const T &dt, const int &step){
    file << dt * step << " ";
    for (T d: distances){
        file << averageNeighbors(particles, d) << " ";
    }
    file << "\n";
}

template<typename T>
void writeAverageNeighborToFile(const std::vector<Particle<T>> &particles, const  T &distance, std::ostream &file, const T &dt, const int &step){
    file << dt * step << " " << averageNeighbors(particles, distance) << "\n";
}