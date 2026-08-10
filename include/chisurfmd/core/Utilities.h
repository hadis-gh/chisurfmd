#pragma once

#include <vector>
#include <cmath>
#include <boost/program_options.hpp>

#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/ParticleDot.h"

namespace po = boost::program_options;

template<typename T>
std::vector<T> linspace(const T& start, const T& end, const int& points) {
    std::vector<T> result;
    T step = (end - start) / (points - 1);
    for (int i = 0; i < points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}

// ================================== order calculation (neighbors + orientation) ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
T calAverageNeighbors(const std::vector<TParticle>& particles, const T& distance) {
    T totalNeighbors = 0;
    for (const auto& p1 : particles) {
        int eachParticleNeighbors = 0;
        for (const auto& p2 : particles) {
            if (p1.position != p2.position && (p1.position - p2.position).abs() <= distance) {
                ++eachParticleNeighbors;
            }
        }
        totalNeighbors += eachParticleNeighbors;
    }
    return totalNeighbors / particles.size();
}

template<typename TParticle, typename T = typename TParticle::value_type>
auto calAveNeighborList(const std::vector<TParticle>& particles, const std::vector<T>& distances, T boxSize) {
    T totalNeighbors = 0;
    std::vector<T> neighborsList;

    for (auto d : distances) {
        totalNeighbors = 0;

        for (const auto& p1 : particles) {
            int eachParticleNeighbors = 0;

            for (const auto& p2 : particles) {
                if (p1.position != p2.position) {
                    // Apply periodic boundary conditions using minimum image convention
                    auto delta = p1.position - p2.position;
                    delta[0] -= boxSize * round(delta[0] / boxSize);
                    delta[1] -= boxSize * round(delta[1] / boxSize);

                    if (delta.abs() <= d) {
                        ++eachParticleNeighbors;
                    }
                }
            }

            totalNeighbors += eachParticleNeighbors;
        }
        
        neighborsList.push_back(totalNeighbors / particles.size());
    }

    return neighborsList;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calOrientationOrder(const std::vector<TParticle>& particles) {

    if constexpr (degreesOfFreedom<TParticle>() <= 2) {
        return T{0};
    }

    T cosSum = 0;
    T sinSum = 0;

    for (const auto& p : particles) {
        cosSum += std::cos(2.0 * p.phi);
        sinSum += std::sin(2.0 * p.phi);
    }

    const T N = static_cast<T>(particles.size());

    cosSum /= N;
    sinSum /= N;

    return std::sqrt(cosSum * cosSum + sinSum * sinSum);
}

// ================================== Center of Mass ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
Vec<T> calCOMposition(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comPos;
    T totalMass = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        comPos += p.position * mass;
        totalMass += mass;
    }
    return comPos / totalMass;
}

template<typename TParticle, typename T = typename TParticle::value_type>
Vec<T> calCOMpositionPBC(const std::vector<TParticle>& particles, 
                         const std::vector<Species<T>>& allSpecies,
                         const T& boxSize) {
    Vec<T> comPos;
    T totalMass = 0;
    
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        comPos += p.position * mass;
        totalMass += mass;
    }
    comPos /= totalMass;

    for (size_t i = 0; i < 2; ++i) {
        comPos[i] = std::fmod(comPos[i], boxSize);
        if (comPos[i] < 0) comPos[i] += boxSize;
    }

    return comPos;
}
// ================================== COM angular velocity ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
T calAngularMomentum2D(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC) {
    T L = 0;
    auto comPos = calCOMpositionPBC(particles, allSpecies, boxPBC);
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        Vec<T> r_com = p.position - comPos;
        L += mass * (r_com[0] * p.velocity[1] - r_com[1] * p.velocity[0]);
    }
    return L;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calMomentOfInertia2D(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const Vec<T>& comPos) {
    T I = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        Vec<T> r_com = p.position - comPos;
        I += mass * r_com.abs2();
    }
    return I;
}

template<typename TParticle, typename T = typename TParticle::value_type>
void removeCOMvelocityRotation2D(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& areaL) {
    Vec<T> comPos = calCOMpositionPBC(particles, allSpecies, areaL);
    T angularMomentum = calAngularMomentum2D(particles, allSpecies, areaL);
    T momentOfInertia = calMomentOfInertia2D(particles, allSpecies, comPos);

    T angularVelocity = angularMomentum / momentOfInertia;

    for (auto& p : particles) {
        Vec<T> r_com = p.position - comPos;
        Vec<T> v_rot = {(-r_com[1] * angularVelocity, r_com[0] * angularVelocity)};
        p.v -= v_rot;
    }
}

template<typename TParticle, typename T = typename TParticle::value_type>
void removeCOMvelocityRotation2D_wholeCenter(std::vector<TParticle>& particles, 
                                            const std::vector<Species<T>>& allSpecies, 
                                            const T& areaL) {
    Vec<T> comPos{{areaL/2, areaL/2}};
    T angularMomentum = calAngularMomentum2D(particles, allSpecies, areaL);
    T momentOfInertia = calMomentOfInertia2D(particles, allSpecies, comPos);

    T angularVelocity = angularMomentum / momentOfInertia;

    for (auto& p : particles) {
        Vec<T> r_com = p.position - comPos;
        Vec<T> v_rot = {(-r_com[1] * angularVelocity, r_com[0] * angularVelocity)};
        p.v -= v_rot;
    }
}

// ================================== COM linear velocity ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
Vec<T> calCOMVelocity(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> totalMomentum = {{0.0, 0.0}};
    T totalMass = 0.0;
    for (const auto& p : particles) {
        T mass = allSpecies[p.species].mass;
        totalMomentum += mass * p.v;
        totalMass += mass;
    }
    return totalMomentum / totalMass;
}

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<T> calCOMvelocity(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> totalMomentum = {{0.0, 0.0}};
    std::vector<T> totalMomentumVec;
    T totalMass = 0.0;
    for (int i = 0; i < 2; ++i) {
        for (const auto& p : particles) {
            T mass = allSpecies[p.species].mass;
            totalMomentum[i] += mass * p.velocity[i];
            totalMass += mass;
        }
        totalMomentumVec.push_back(totalMomentum[i] / totalMass);
    }

    return totalMomentumVec;
}

// ================================== remove COM angular velocity ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
void removeCOMVelocity(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comVel = calCOMVelocity(particles, allSpecies);
    for (auto& p : particles) {
        p.v -= comVel;
    }
}

// ================================== Move particles to center of box ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
void moveParticlesToCenter(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& L)
{
    Vec<T> comPos = calCOMpositionPBC(particles, allSpecies, L);
    Vec<T> shift = -comPos + L/2;
    
    for (auto &p: particles){
        p.position += shift;
        for (size_t i = 0; i < 2; ++i) {
            p.position[i] = std::fmod(p.position[i], L);
            if (p.position[i] < 0) p.position[i] += L;
        }
    }     
}

// ================================== Fix particles around the moving one ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
void applyFixRadius(std::vector<TParticle>& particles, const T& fixRadius, const T& areaL)
{
    Vec<T> boxCenter {{areaL / 2, areaL / 2}};
    for (auto& p : particles) {
        if ((p.position - boxCenter).abs2() > fixRadius * fixRadius) {
            p.fixed = true;
        }
    }
}