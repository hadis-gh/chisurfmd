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
auto calAveNeighborList(const std::vector<TParticle>& particles, const std::vector<T>& distances) {
    T totalNeighbors = 0;
    std::vector<T> neighborsList;
    for (auto &d : distances) {
        T totalNeighbors = 0;
        for (const auto& p1 : particles) {
            int eachParticleNeighbors = 0;
            for (const auto& p2 : particles) {
                if (p1.r != p2.r && (p1.r - p2.r).abs() <= d) {
                    ++eachParticleNeighbors;
                }
            }
            totalNeighbors += eachParticleNeighbors;
        }
        neighborsList.push_back(totalNeighbors / particles.size());
    }
    return neighborsList;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calPositionalOrder (const std::vector<TParticle>& particles, const T& neighborCutoff) {
    //implement later
    return neighborCutoff;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calOrientationalOrder(const std::vector<TParticle>& particles) {
    T orientationalOrder;
    size_t count = 0;


    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            auto pr1 = getGeneralizedPositions(particles[i]);
            auto pr2 = getGeneralizedPositions(particles[j]);

            if (pr1.size() > 2) {
                T deltaPhi = pr1[2] - pr2[2];
                deltaPhi = std::fmod(deltaPhi + 2*M_PI, 2*M_PI);
                if (deltaPhi > M_PI) deltaPhi -= 2*M_PI;
    
                orientationalOrder += std::cos(2.0 * deltaPhi * M_PI / M_PI);
                count++;
            } else {
                orientationalOrder = 0;
                count = 1;
            }
        }
    }
    return orientationalOrder/count;
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
Vec<T> calCOMpositionPBC(const std::vector<TParticle>& particles, 
                         const std::vector<Species<T>>& allSpecies,
                         const T& boxSize) {
    Vec<T> comPos;
    T totalMass = 0;
    Vec<T> prevPos = particles[0].r;
    
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        Vec<T> adjustedPos = p.r;
        for (size_t i = 0; i < 2; ++i) {
            while (adjustedPos[i] - prevPos[i] > boxSize/2) adjustedPos[i] -= boxSize;
            while (adjustedPos[i] - prevPos[i] < -boxSize/2) adjustedPos[i] += boxSize;
        }
        comPos += adjustedPos * mass;
        totalMass += mass;
        prevPos = adjustedPos;
    }
    comPos /= totalMass;
    
    for (size_t i = 0; i < 2; ++i) {
        comPos[i] = std::fmod(comPos[i], boxSize);
        if (comPos[i] < 0) comPos[i] += boxSize;
    }
    
    return comPos;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calAngularMomentum2D(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC) {
    T L = 0;
    auto comPos = calCOMpositionPBC(particles, allSpecies, boxPBC);
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
void removeCOMvelocityRotation2D(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& areaL) {
    Vec<T> comPos = calCOMpositionPBC(particles, allSpecies, areaL);
    T angularMomentum = calAngularMomentum2D(particles, allSpecies, areaL);
    T momentOfInertia = calMomentOfInertia2D(particles, allSpecies, comPos);

    T angularVelocity = angularMomentum / momentOfInertia;

    for (auto& p : particles) {
        Vec<T> r_com = p.r - comPos;
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
        Vec<T> r_com = p.r - comPos;
        Vec<T> v_rot = {(-r_com[1] * angularVelocity, r_com[0] * angularVelocity)};
        p.v -= v_rot;
    }
}

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
            totalMomentum[i] += mass * p.v[i];
            totalMass += mass;
        }
        totalMomentumVec.push_back(totalMomentum[i] / totalMass);
    }

    return totalMomentumVec;
}

template<typename TParticle, typename T = typename TParticle::value_type>
void removeCOMVelocity(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comVel = calCOMVelocity(particles, allSpecies);
    for (auto& p : particles) {
        p.v -= comVel;
    }
}