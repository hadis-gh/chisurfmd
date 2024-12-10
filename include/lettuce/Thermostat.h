#pragma once

#include <vector>
#include <cmath>
#include <string>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/constants.h"
#include "lettuce/LennardJones.h"

template<typename TParticle, typename T = typename TParticle::value_type>
class AndersenThermostat {
public:
    AndersenThermostat(T dt, T collisionFrequency, T desiredTemperature, std::mt19937 gen)
        : dt(dt), collisionFrequency(collisionFrequency), desiredTemperature(desiredTemperature),
          dist(0.0, 1.0), maxwellDist(0.0, std::sqrt(constants::boltzmann * desiredTemperature)) {}

    void operator () (std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
        for (auto& p : particles) {
            if (dist(gen) < collisionFrequency * dt) {
                const T mass = allSpecies[p.species].mass;
                const T MI = allSpecies[p.species].momentOfInertia;
                auto generalizedVel = getGeneralizedVelocities(p);

                for (size_t i = 0; i < generalizedVel.size(); ++i) {
                    if (i < 2) {
                        generalizedVel[i] = maxwellDist(gen) / std::sqrt(mass);
                    } else {
                        generalizedVel[i] = maxwellDist(gen) / std::sqrt(MI);
                    }
                }

                setGeneralizedVelocities(p, generalizedVel);
            }
        }
    }

private:
    T dt;
    T collisionFrequency;
    T desiredTemperature;
    std::mt19937 gen;
    std::uniform_real_distribution<T> dist;
    std::normal_distribution<T> maxwellDist;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class BerendsenThermostat {
public:
    BerendsenThermostat(T dt, T desiredTemperature, T relaxationTime)
        : dt(dt), desiredTemperature(desiredTemperature), relaxationTime(relaxationTime) {}

    void operator () (std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) const {
        const T currentTemperature = calInternalTemperature(particles, allSpecies); 
        const auto lambda = std::sqrt(1 + dt / relaxationTime * (desiredTemperature / currentTemperature - 1));
        rescaleVelocity(particles, allSpecies, lambda);
    }

private:
    T dt;
    T relaxationTime;
    T desiredTemperature;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class VelocityScalingThermostat {
public:
    VelocityScalingThermostat(T desiredTemperature)
        : desiredTemperature(desiredTemperature) {}

    void operator () (std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) const {
        const T currentTemperature = calInternalTemperature(particles, allSpecies); 
        const auto lambda = std::sqrt(desiredTemperature / currentTemperature);
        rescaleVelocities(particles, allSpecies, lambda);
    }

private:
    T desiredTemperature;
};

template<typename TParticle, typename T = typename TParticle::value_type>
void rescaleVelocities(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T lambda) {
    for (auto& p : particles) {
        auto generalizedVel = getGeneralizedVelocities(p);
        const T mass = allSpecies[p.species].mass;
        const T MI = allSpecies[p.species].momentOfInertia;

        for (size_t i = 0; i < generalizedVel.size(); ++i) {
            if (i < 2) {
                generalizedVel[i] *= std::sqrt(lambda / mass);
            } else {
                generalizedVel[i] *= std::sqrt(lambda / MI);
            }
        }
        setGeneralizedVelocities(p, generalizedVel);
    }
}
