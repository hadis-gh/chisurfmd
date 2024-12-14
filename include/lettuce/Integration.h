#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

#include <adios2.h>

#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/LennardJones.h"
#include "lettuce/LennardJonesOriented.h"
#include "lettuce/Thermostat.h"

template<typename TParticle, typename Force, typename T = typename TParticle::value_type>
void EulerSymplecticStep(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto r = getGeneralizedPositions(particles[i]);
        const auto v = getGeneralizedVelocities(particles[i]);
        setGeneralizedVelocities(particles[i], v + accelerations[i] * dt);
        setGeneralizedPositions(particles[i], r + particles[i].v * dt);
        implementPBC(particles[i], boxPBC);
    }
}

template<typename TParticle, typename Force, typename T = typename TParticle::value_type>
void EulerStep(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto r = getGeneralizedPositions(particles[i]);
        const auto v = getGeneralizedVelocities(particles[i]);
        setGeneralizedPositions(particles[i], r + v * dt);
        setGeneralizedVelocities(particles[i], v + accelerations[i] * dt);
        implementPBC(particles[i], boxPBC);
    }
}

template<typename TParticle, typename Force, typename T = typename TParticle::value_type>
void VelocityVerletStep(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto old_accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto r = getGeneralizedPositions(particles[i]);
        const auto v = getGeneralizedVelocities(particles[i]);

        setGeneralizedPositions(particles[i], r + v * dt + old_accelerations[i] * dt * dt / 2.0);
        implementPBC(particles[i], boxPBC);
    }

    const auto new_accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto v = getGeneralizedVelocities(particles[i]);
        setGeneralizedVelocities(particles[i], v + (old_accelerations[i] + new_accelerations[i]) * dt / 2.0);
    }
}

template<typename TParticle, typename Integrator, typename Force, typename T = typename TParticle::value_type>
void integrate(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, T dt, T Time, const T& boxPBC, Force&& force, Integrator&& integrator) {
    int numSteps = static_cast<int>(Time / dt);

    for (int i = 0; i < numSteps; ++i) {
        integrator(particles, allSpecies, dt, boxPBC, std::forward<Force>(force));
    }
}

// template<typename TParticle>
// void capVelocity(TParticle& particle, const T max_velocity) {
//     auto v = getGeneralizedVelocities(particle);
//     T speed = v.abs();
//     if (speed > max_velocity) {
//         v *= max_velocity / speed;
//         setGeneralizedVelocities(particle, v);
//     }
// }
