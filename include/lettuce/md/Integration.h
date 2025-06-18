#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>
#include <cstdlib>

#include <adios2.h>

#include "lettuce/core/Vec.h"
#include "lettuce/core/Circle.h"
#include "lettuce/md/Thermostat.h"

template<typename TParticle, typename Force, typename T = typename TParticle::value_type>
void EulerSymplecticStep(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto r = getGeneralizedPositions(particles[i]);
        const auto v = getGeneralizedVelocities(particles[i]);
        setGeneralizedVelocities(particles[i], v + accelerations[i] * dt);
        setGeneralizedPositions(particles[i], r + v * dt);
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

template<typename T, typename Particle>
void capVelocity(std::vector<Particle>& particles, const std::vector<T>& capV) {
    for (auto& particle : particles) {
        auto velocities = getGeneralizedVelocities(particle);

        const T translationalCap = capV[0];
        for (size_t i = 0; i < 2; ++i) {
            if (std::abs(velocities[i]) > translationalCap) {
                velocities[i] = (velocities[i] > 0 ? 1 : -1) * translationalCap;
            }
        }

        if (velocities.size() > 2) {
            const T rotationalCap = capV[1];
            if (std::abs(velocities[2]) > rotationalCap) {
                velocities[2] = (velocities[2] > 0 ? 1 : -1) * rotationalCap;
            }
        }

        setGeneralizedVelocities(particle, velocities);
    }
}

