#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/ParticleDot.h"

template<typename T>
std::vector<T> linspace(const T& start, const T& end, const int& points) {
    std::vector<T> result;
    T step = (end - start) / (points - 1);
    for (int i = 0; i < points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calAverageNeighbors(const std::vector<TParticle>& particles, const T& distance) {
    T totalNeighbors = 0;
    for (const auto& p1 : particles) {
        int eachParticleNeighbors = 0;
        for (const auto& p2 : particles) {
            if (p1.r != p2.r && (p1.r - p2.r).abs() <= distance) {
                ++eachParticleNeighbors;
            }
        }
        totalNeighbors += eachParticleNeighbors;
    }
    return totalNeighbors / particles.size();
}

template<typename TParticle, typename T = typename TParticle::value_type>
Vec<T> calCOMposition(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comPos;
    T totalMass = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        comPos += p.r * mass;
        totalMass += mass;
    }
    return comPos / totalMass;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calAngularMomentum2D(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const Vec<T>& comPos) {
    T L = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        Vec<T> r_com = p.r - comPos;
        L += mass * (r_com[0] * p.v[1] - r_com[1] * p.v[0]);
    }
    return L;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calMomentOfInertia2D(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const Vec<T>& comPos) {
    T I = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        Vec<T> r_com = p.r - comPos;
        I += mass * r_com.abs2();
    }
    return I;
}

template<typename TParticle, typename T = typename TParticle::value_type>
void removeCOMvelocityRotation2D(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comPos = calCOMposition(particles, allSpecies);
    T angularMomentum = calAngularMomentum2D(particles, allSpecies, comPos);
    T momentOfInertia = calMomentOfInertia2D(particles, allSpecies, comPos);

    T angularVelocity = angularMomentum / momentOfInertia;
    // there is possiblity to remove mass from these eqs - try to simplify it

    for (auto& p : particles) {
        Vec<T> r_com = p.r - comPos;
        Vec<T> v_rot = {(-r_com[1] * angularVelocity, r_com[0] * angularVelocity)};
        p.v -= v_rot;
    }
}

template<typename TParticle, typename T = typename TParticle::value_type>
Vec<T> calCOMVelocity(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> totalMomentum;
    T totalMass = 0.0;
    for (const auto& p : particles) {
        T mass = allSpecies[p.species].mass;
        totalMomentum += mass * p.v;
        totalMass += mass;
    }
    return totalMomentum / totalMass;
}

template<typename TParticle, typename T = typename TParticle::value_type>
void removeCOMVelocity(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comVel = calCOMVelocity(particles, allSpecies);
    for (auto& p : particles) {
        p.v -= comVel;
    }
}