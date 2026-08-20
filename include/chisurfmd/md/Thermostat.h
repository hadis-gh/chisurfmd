#pragma once

#include <vector>
#include <cmath>
#include <string>
#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/Circle.h"
#include "chisurfmd/constants.h"

template<typename TParticle, typename T = typename TParticle::value_type>
class AndersenThermostat {
public:
    AndersenThermostat(T dt, T collisionFrequency, T desiredTemperature, std::mt19937 gen)
        : dt(dt), collisionFrequency(collisionFrequency), desiredTemperature(desiredTemperature),
          dist(0.0, 1.0),
          translationalDist(0.0, std::sqrt(constants::boltzmann * desiredTemperature)),
          rotationalDist(0.0, std::sqrt(constants::boltzmann * desiredTemperature)),
          gen(gen) {}

    void operator () (std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
        for (auto& p : particles) {
            if (dist(gen) < collisionFrequency * dt) {
                const T mass = allSpecies[p.species].mass;
                const T MI = allSpecies[p.species].momentOfInertia;
                auto generalizedVel = getGeneralizedVelocities(p);

                for (size_t i = 0; i < 2; ++i) {
                    generalizedVel[i] = translationalDist(gen) * std::sqrt(1.0 / mass);
                }

                if (generalizedVel.size() > 2) {
                    generalizedVel[2] = rotationalDist(gen) * std::sqrt(1.0 / MI);  //test to remove this part and see if it works
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
    std::normal_distribution<T> translationalDist;
    std::normal_distribution<T> rotationalDist;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class BerendsenThermostat {
public:
    BerendsenThermostat(T dt, T desiredTemperature, T relaxationTime)
        : dt(dt), desiredTemperature(desiredTemperature), relaxationTime(relaxationTime) {}

    void operator () (std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) const {
        const T currentTemperature = calInternalTemperature(particles, allSpecies); 
        const auto lambda = std::sqrt(1 + dt / relaxationTime * (desiredTemperature / currentTemperature - 1));
        rescaleVelocities(particles, allSpecies, lambda);
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

// ================================== Reset Velocities after DLA deposition ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
void resetVelocitiesRandom(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, T desiredTemperature, std::mt19937& gen) {
    
    std::normal_distribution<T> translationalDist(0.0, std::sqrt(constants::boltzmann * desiredTemperature));
    std::normal_distribution<T> rotationalDist(0.0, std::sqrt(constants::boltzmann * desiredTemperature));

    for (auto& p : particles) {
        auto generalizedVel = getGeneralizedVelocities(p);
        const T mass = allSpecies[p.species].mass;
        const T MI = allSpecies[p.species].momentOfInertia;
        
        for (size_t i = 0; i < generalizedVel.size(); ++i) {
            if (i < 2) {
                generalizedVel[i] = translationalDist(gen) * std::sqrt(1.0 / mass);
            } else {
                generalizedVel[2] = rotationalDist(gen) * std::sqrt(1.0 / MI);
            }
        }
        setGeneralizedVelocities(p, generalizedVel);
    }
}