#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/LennardJones.h"

enum IntegrationMethod {
    EULER,
    VELOCITY_VERLET
};

template<typename T, typename Force>        //Euler Integration
void updateStateEuler(std::vector<Particle<T>> &particles, const T &dt, Force &&force){
    const auto accelerations = calAllAccelerations(particles, std::forward<Force>(force));
    for (unsigned int a = 0; a < particles.size(); ++a) {
        auto& p = particles[a];
        p.v += accelerations[a]*dt;
        p.r += p.v * dt;
    }
}

template<typename T, typename Force>        //Velocity Vernel Integration
void updateStateVV(std::vector<Particle<T>> &particles, const T &dt, Force &&force) {
    const auto old_accelerations = calAllAccelerations(particles, std::forward<Force>(force));
    for (unsigned int a = 0; a < particles.size(); ++a) {
        auto& p = particles[a];
        p.r += p.v * dt + old_accelerations[a]/2 * dt * dt;
    }
    const auto accelerations = calAllAccelerations(particles, std::forward<Force>(force));
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].v += (old_accelerations[i] + accelerations[i])/2 * dt;
    }
}

template<typename T>
void writeToFile(const std::vector<Particle<T>> &particles, std::ostream &file) {
    for (const auto &p : particles) {
        file << p.r[0] << " " << p.r[1] << " ";
    }
    file << "\n";
}

template<typename T>
void integrate(std::vector<Particle<T>>& particles, T dt, T Time, LennardJonesForce<T>& LJForce, IntegrationMethod method) {
    std::string filename;
    if (method == EULER) {
        filename = "particlesPosMD_Euler.dat";
    } else if (method == VELOCITY_VERLET) {
        filename = "particlesPosMD_VV.dat";
    }
    std::ofstream file(filename);
    writeToFile(particles, file);

    int numSteps = static_cast<int>(Time / dt);
    if (method == EULER) {
        for (int i = 0; i < numSteps; ++i) {
            updateStateEuler(particles, dt, LJForce);
            writeToFile(particles, file);
        }
    } else if (method == VELOCITY_VERLET) {
        for (int i = 0; i < numSteps; ++i) {
            updateStateVV(particles, dt, LJForce);
            writeToFile(particles, file);
        }
    }
    file.close();
}