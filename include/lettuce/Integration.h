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
void updateStateEuler(std::vector<Particle<T>> &particles, std::vector<Species<T>> allSpecies, const T &dt, Force &&force){
    const auto accelerations = calAllAccelerations(particles, allSpecies, std::forward<Force>(force));
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        p.v += accelerations[i]*dt;
        p.r += p.v * dt;
    }
}

template<typename T, typename Force>        //Velocity Vernel Integration
void updateStateVV(std::vector<Particle<T>> &particles, std::vector<Species<T>>& allSpecies, const T &dt, Force &&force) {
    const auto old_accelerations = calAllAccelerations(particles, allSpecies, std::forward<Force>(force));
    for (size_t i = 0; i < particles.size(); ++i) {
        auto& p = particles[i];
        p.r += p.v * dt + old_accelerations[i]/2 * dt * dt;
    }
    const auto accelerations = calAllAccelerations(particles, allSpecies, std::forward<Force>(force));
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].v += (old_accelerations[i] + accelerations[i])/2 * dt;
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
void writeKineticEToFile(const std::vector<Particle<T>> &particles, std::vector<Species<T>>& allSpecies, std::ostream &file) {
    for (const auto &p : particles) {
        file << 0.5 * allSpecies[p.species].mass * p.v.abs2() * p.v.abs2() << " ";
    }
    file << "\n";
}

// template<typename T, typename Force>
// class integrate{
// public:
//     integrate(std::vector<Particle<T>>& particles, T dt, T Time, LennardJonesForce<T>& LJForce, IntegrationMethod method)

// private:

// }

template<typename T>
void integrate(std::vector<Particle<T>>& particles, std::vector<Species<T>>& allSpecies, T dt, T Time, LennardJonesForce<T>& LJForce, IntegrationMethod method) {
    std::string filenamePos;
    std::string filenameKineticE;

    if (method == EULER) {
        filenamePos = "particlesPosMD_Euler.dat";
        filenameKineticE = "particlesKineticE_Euler.dat";
    } else if (method == VELOCITY_VERLET) {
        filenamePos = "particlesPosMD_VV.dat";
        filenameKineticE = "particlesKineticE_VV.dat";
    }
    std::ofstream filePos(filenamePos);
    std::ofstream fileKE(filenameKineticE);

    writePositionToFile(particles, filePos);
    writeKineticEToFile(particles, allSpecies, fileKE);

    int numSteps = static_cast<int>(Time / dt);
    if (method == EULER) {
        for (size_t i = 0; i < numSteps; ++i) {
            updateStateEuler(particles, allSpecies, dt, LJForce);
            writePositionToFile(particles, filePos);
            writeKineticEToFile(particles, allSpecies, fileKE);
        }
    } else if (method == VELOCITY_VERLET) {
        for (size_t i = 0; i < numSteps; ++i) {
            updateStateVV(particles, allSpecies, dt, LJForce);
            writePositionToFile(particles, filePos);
            writeKineticEToFile(particles, allSpecies, fileKE);
        }
    }
    filePos.close();
    fileKE.close();
}