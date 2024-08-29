#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"


template<typename T>
struct Species
{
    T mass;
    T radius;
    T I = 0.5 * mass * radius * radius;
};

template<typename T>
struct Particle
{
    Vec<T> r, v;
    T w;
    unsigned int species;
};

template<typename T, typename Force>
Vec<T> calForceTwo(const Particle<T> &p1, const Particle<T> &p2, const T& boxPBC, Force &&force){
    Vec<T> dr = p2.r - p1.r;
    
    for (int i = 0; i < 2; ++i) {
        if (dr[i] > boxPBC / 2) { dr[i] -= boxPBC; }
        else if (dr[i] < -boxPBC / 2) { dr[i] += boxPBC; }
    }
    
    T r = dr.abs();
    if (r == 0) return {{0, 0}};
    
    T f = force(r);
    return f * dr / r;
}

template<typename T, typename Force>
Vec<T> calTotalForce(const Particle<T> &p1, const std::vector<Particle<T>> &particles, const T& boxPBC, Force &&force){
    Vec<T> f;
    for (auto &p : particles){
        if (p1.r != p.r){
            f += calForceTwo(p, p1, boxPBC, std::forward<Force>(force));
        }
    }
    return f;
}

template<typename T, typename Force>
Vec<T> calAccelaration(const Particle<T> &p1, const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Force &&force){
    T mass = allSpecies[p1.species].mass;
    return calTotalForce(p1, particles, boxPBC, std::forward<Force>(force))/ mass;
}

template<typename T, typename Force>
std::vector<Vec<T>> calAllAccelerations(const std::vector<Particle<T>> &particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Force &&force) {
    std::vector<Vec<T>> accelerations(particles.size());
    for (unsigned int i = 0; i < particles.size(); ++i) {
        accelerations[i] = calAccelaration(particles[i], particles, allSpecies, boxPBC, std::forward<Force>(force));
    }
    return accelerations;
}

template<typename T>
T calculateTorque(const Particle<T>& p, const Vec<T>& force) {
    return p.r[0] * force[1] - p.r[1] * force[0];
}
