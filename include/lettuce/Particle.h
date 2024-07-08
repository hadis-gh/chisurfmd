#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"

const double cutoff = 200;

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
        if ((p1.r - p.r) <= cutoff){
            if (p1.r != p.r){
                f += calForceTwo(p, p1, force);
            }
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
    for (unsigned int a = 0; a < particles.size(); ++a) {
        accelerations[a] = calAccelaration(particles[a], particles, force);
    }
    return accelerations;
}