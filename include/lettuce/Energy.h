#pragma once

#include <vector>
#include <cmath>
#include <iomanip>

#include "lettuce/Vec.h"
#include "lettuce/constants.h"

// template<typename TParticle, typename T = typename TParticle::value_type>
// T calInternalTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
//     T internalKE = calInternalKineticEnergy(particles, allSpecies);
//     return (2 * internalKE) / (constants::boltzmann * particles.size() * 3.0);
// }

template<typename TParticle, typename T = typename TParticle::value_type>
T calInternalTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T translationalKE = 0.0;
    T rotationalKE = 0.0;
    size_t numParticles = particles.size();
    
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        const T MI = allSpecies[p.species].momentOfInertia;
        auto vel = getGeneralizedVelocities(p);

        for (size_t i = 0; i < 2; ++i) {
            translationalKE += 0.5 * mass * vel[i] * vel[i];
        }

        if (vel.size() > 2) {
            rotationalKE += 0.5 * MI * vel[2] * vel[2];
        }
    }

    // Each degree of freedom contributes 1/2 kB T to the energy: equipartition theorem
    T totalDOF = 2 * numParticles;  // 2 translational DOF per particle
    if (DegreesOfFreedom<TParticle>::degreesOfFreedom() > 2) {
        totalDOF += numParticles;    // Add 1 rotational DOF per particle if present
    }

    return 2.0 * (translationalKE + rotationalKE) / (constants::boltzmann * totalDOF);
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calRawTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T rawKE = calRawKineticEnergy(particles, allSpecies);
    return (2 * rawKE) / (constants::boltzmann * particles.size() * 3.0);
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calInternalKineticEnergy(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comVel = calCOMVelocity(particles, allSpecies);
    T totalKE = 0.0;
    for (const auto& p : particles) {
        Vec<T> relativeVelocity = p.v - comVel;
        totalKE += calParticleRelativeKineticEnergy(p, allSpecies, relativeVelocity);
    }
    return totalKE;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calRawKineticEnergy(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T totalKE = 0.0;
    for (const auto& p : particles) {
        totalKE += calParticleKineticEnergy(p, allSpecies);
    }
    return totalKE;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calParticleKineticEnergy(const TParticle& p, const std::vector<Species<T>>& allSpecies) {
    T kineticEnergy = 0.0;
    const T mass = allSpecies[p.species].mass;
    const T MI = allSpecies[p.species].momentOfInertia;
    auto vel = getGeneralizedVelocities(p);

    for (size_t i = 0; i < vel.size(); ++i) {
        if (i < 2) {
            kineticEnergy += 0.5 * mass * vel[i] * vel[i];
        } else {
            kineticEnergy += 0.5 * MI * vel[i] * vel[i];
        }
    }
    return kineticEnergy;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calParticleRelativeKineticEnergy(const TParticle& p, const std::vector<Species<T>>& allSpecies, const Vec<T>& relativeVelocity) {
    return 0.5 * allSpecies[p.species].mass * relativeVelocity.abs2();
}

template<typename TParticle, typename Potential, typename T = typename TParticle::value_type>
T calPotentialEnergy(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Potential&& potential) {
    std::vector<T> ePot(particles.size());
    for (auto &p1 : particles) {
        double pot = 0.;
        for (auto &p : particles) {
            if (p1.r != p.r) {
                Vec<T> dr = p.r - p1.r;

                for (int i = 0; i < 2; ++i) {
                    if (dr[i] > boxPBC / 2) { dr[i] -= boxPBC; }
                    else if (dr[i] < -boxPBC / 2) { dr[i] += boxPBC; }
                }

                const T r = dr.abs();
                if (r == 0) continue;

                if constexpr (std::is_same_v<TParticle, ParticleOriented<T>>) {
                    T deltaPhi = p.phi - p1.phi;
                    pot += potential(r, deltaPhi);
                } else {
                    pot += potential(r);
                }
            }
        }
        ePot.push_back(pot);
    }
    return std::accumulate(ePot.begin(), ePot.end(), 0.0);
}

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<T> calKineticEnergy(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    T totalKE_x = 0.0;
    T totalKE_y = 0.0;
    T totalKE_phi = 0.0;

    size_t numParticles = particles.size();

    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        const T MOI = allSpecies[p.species].momentOfInertia;

        auto vel = getGeneralizedVelocities(p);

        totalKE_x += 0.5 * mass * vel[0] * vel[0];
        totalKE_y += 0.5 * mass * vel[1] * vel[1];

        if (vel.size() > 2) {
            totalKE_phi += 0.5 * MOI * vel[2] * vel[2];
        }
    }

    T avgKE_x = totalKE_x / numParticles;
    T avgKE_y = totalKE_y / numParticles;
    T avgKE_phi = (totalKE_phi > 0) ? (totalKE_phi / numParticles) : 0.0;

    return {{avgKE_x, avgKE_y, avgKE_phi}};
}