#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"


template<typename T>
struct Particle{
    Vec<T> r;
    Vec<T> v;
    T mass;
};

template<typename T, typename Force>
Vec<T> calForceTwo(const Particle<T> &p1, const Particle<T> &p2, Force force){
    T r = (p2.r - p1.r).abs2();
    if (r == 0) return {{0, 0}};
    
    T f = force(r);
    return f * (p2.r - p1.r)/r;
}

template<typename T, typename Force>
Vec<T> calTotalForce(const Particle<T> &p1, const std::vector<Particle<T>> &particles, Force force){
    Vec<T> f;
    for (auto &p : particles){
        if (p1.r != p.r){
            f += calForceTwo(p, p1, force);
        }
    }
    return f;
}

template<typename T, typename Force>
Vec<T> calAccelaration(const Particle<T> &p1, const std::vector<Particle<T>> &particles, Force force){
    return calTotalForce(p1, particles, force)/ p1.mass;
}

template<typename T, typename Force>
std::vector<Vec<T>> calAllAccelerations(std::vector<Particle<T>> &particles, Force force) {
    std::vector<Vec<T>> accelerations(particles.size());
    for (unsigned int i = 0; i < particles.size(); ++i) {
        accelerations[i] = calAccelaration(particles[i], particles, force);
    }
    return accelerations;
}

template<typename T>
std::vector<T> calculateKineticEnergy (std::vector<Particle<T>> &particles, const T &Time, const T &dt){
    const T steps = Time/ dt;
    std::vector<T> kineticEnergy(steps);
    
    for (int i = 0; i < steps; ++i){
        for (int j= 0; j < particles.size(); ++j){
            kineticEnergy[i] = 0.5 * particles[0].mass * particles[j].v * particles[j].v;
        }
    }
    return kineticEnergy;
}

template<typename T, typename Force>
std::vector<T> calculatePotentialEnergy (std::vector<Particle<T>> &particles, const T &Time, const T &dt, Force potential){
    const T steps = Time/ dt;
    std::vector<Vec<T>> accelarations(steps);
    std::vector<T> potentialEnergy(steps);

    accelarations = calAllAccelerations(particles, potential);
    for (int i = 0; i < steps; ++i){
        potentialEnergy[i] = accelarations[i].abs2();
    }
    return potentialEnergy;
}

template<typename T, typename Force>
std::vector<T> calculateTotalEnergy (std::vector<Particle<T>> &particles, const T &Time, const T &dt, Force potentialType){
    const T steps = Time/ dt;    
    std::vector<T> potentialEnergy(steps);
    std::vector<T> kineticEnergy(steps);
    std::vector<T> totalEnergy(steps);

    potentialEnergy = calculatePotentialEnergy (particles, Time, dt, potentialType);
    kineticEnergy = calculateKineticEnergy (particles, Time, dt);

    for (int i = 0; i < steps; ++i){
        totalEnergy[i] =  kineticEnergy[i] + potentialEnergy[i];
    }

    return totalEnergy;
}