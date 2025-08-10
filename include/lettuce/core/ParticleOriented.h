#pragma once

#include <vector>
#include <cmath>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"

template<typename T>
struct ParticleOriented : public ParticleDot<T>
{
    using typename ParticleDot<T>::value_type;
    using ParticleDot<T>::ParticleDot;

    value_type phi, omega=0;

    ParticleOriented() : ParticleDot<T>(), phi(0), omega(0) {}

    ParticleOriented(value_type phi, value_type omega)
        : ParticleDot<T>(), phi(phi), omega(omega) {}

    ParticleOriented(unsigned int speciesIndex, value_type phi, value_type omega)
        : ParticleDot<T>(speciesIndex), phi(phi), omega(omega) {}
};

template<typename T>
struct CreateRandomParticle<ParticleOriented<T>>
{
    static ParticleOriented<T> createRandomParticle(std::mt19937& gen)
    {

        ParticleDot<T> baseParticle = CreateRandomParticle<ParticleDot<T>>::createRandomParticle(gen);

        ParticleOriented<T> p;
        p.r[0] = baseParticle.r[0];
        p.r[1] = baseParticle.r[1];

        std::uniform_real_distribution<T> randomPos(0, 1);
        // p.phi = M_PI / 2;
        p.phi = randomPos(gen) * 2 * M_PI - M_PI;

        return p;
    }
};

template<typename T>
struct CreateTwoParticle<ParticleOriented<T>> {
    static std::vector<ParticleOriented<T>> createTwoParticle(std::mt19937& gen, const T& areaL) {
        auto baseParticles = CreateTwoParticle<ParticleDot<T>>::createTwoParticle(gen, areaL);
        
        std::vector<ParticleOriented<T>> orientedParticles;
        for (auto& baseParticle : baseParticles) {
            ParticleOriented<T> p;
            p.r[0] = baseParticle.r[0];
            p.r[1] = baseParticle.r[1];

            std::uniform_real_distribution<T> randomPos(0, 1);
            p.phi = randomPos(gen) * 2 * M_PI - M_PI;       ///equilibrium phi

            orientedParticles.push_back(p);
        }
        return orientedParticles;
    }
};

template<typename T>
struct DegreesOfFreedom<ParticleOriented<T>>
{
    static constexpr int degreesOfFreedom()
    {
        return 3;
    }
};

// getting generalized positions
template<typename T>
struct GetGeneralizedPositions<ParticleOriented<T>>
{
    static Vec<T,3> getGeneralizedPositions(const ParticleOriented<T>& p)
    {
        Vec<T,3> r;
        r[0] = p.r[0];
        r[1] = p.r[1];
        r[2] = p.phi;
        return r;
    }
};

// setting generalized positions
template<typename T>
struct SetGeneralizedPositions<ParticleOriented<T>>
{
    static void setGeneralizedPositions(ParticleOriented<T>& p, const Vec<T, 3>& new_positions)
    {
        p.r[0] = new_positions[0];
        p.r[1] = new_positions[1];
        p.phi = new_positions[2];
    }
};

// getting generalized velocities
template<typename T>
struct GetGeneralizedVelocities<ParticleOriented<T>>
{
    static Vec<T, 3> getGeneralizedVelocities(const ParticleOriented<T>& p)
    {
        Vec<T, 3> v;
        v[0] = p.v[0];
        v[1] = p.v[1];
        v[2] = p.omega;
        return v;
    }
};

// setting generalized velocities
template<typename T>
struct SetGeneralizedVelocities<ParticleOriented<T>>
{
    static void setGeneralizedVelocities(ParticleOriented<T>& p, const Vec<T, 3>& new_velocities)
    {
        p.v[0] = new_velocities[0];
        p.v[1] = new_velocities[1];
        p.omega = new_velocities[2];
    }
};

// Force calculations for ParticleOriented
template<typename T, typename Force>
Vec<T, 3> calForceTwo(const ParticleOriented<T>& p1, const ParticleOriented<T>& p2, 
                      const T& boxPBC, Force&& force) {
    Vec<T, 2> dr = p2.r - p1.r;
    
    for (int i = 0; i < 2; ++i) {
        if (dr[i] > boxPBC / 2) dr[i] -= boxPBC;
        else if (dr[i] < -boxPBC / 2) dr[i] += boxPBC;
    }
    
    const T r = dr.abs();
    if (r <= std::numeric_limits<T>::epsilon()) return {{0, 0, 0}};
    
    return force(p1, p2, dr, r);
}

template<typename T, typename Force>
Vec<T, 3> calTotalForce(const ParticleOriented<T>& p1, 
                        const std::vector<ParticleOriented<T>>& particles, 
                        const T& boxPBC, Force&& force) {
    Vec<T, 3> totalForce {{0, 0, 0}};
    for (const auto& p : particles) {
        if (&p1 != &p) {
            totalForce += calForceTwo(p, p1, boxPBC, std::forward<Force>(force));
        }
    }
    return totalForce;
}

template<typename T, typename Force>
Vec<T, 3> calAcceleration(const ParticleOriented<T>& p1, 
                          const std::vector<ParticleOriented<T>>& particles, 
                          const std::vector<Species<T>>& allSpecies, 
                          const T& boxPBC, Force&& force) {
    const T mass = allSpecies[p1.species].mass;
    const T MI = allSpecies[p1.species].momentOfInertia;
    Vec<T, 3> totalForce = calTotalForce(p1, particles, boxPBC, std::forward<Force>(force));
    
    return {{
        totalForce[0] / mass,
        totalForce[1] / mass,
        totalForce[2] / MI
    }};
}

template<typename T, typename Force>
std::vector<Vec<T, 3>> calAllAccelerations(const std::vector<ParticleOriented<T>>& particles, 
                                          const std::vector<Species<T>>& allSpecies, 
                                          const T& boxPBC, Force&& force) {
    std::vector<Vec<T, 3>> accelerations;
    accelerations.reserve(particles.size());
    
    for (const auto& p : particles) {
        accelerations.push_back(
            calAcceleration(p, particles, allSpecies, boxPBC, std::forward<Force>(force)));
    }
    return accelerations;
}

template<typename T>
void implementPBC(ParticleOriented<T>& p, const T& boxPBC) {
    implementPBC(static_cast<ParticleDot<T>&>(p), boxPBC);
    
    // Normalize angle to [-π, π)
    p.phi = std::fmod(p.phi + M_PI, 2 * M_PI);
    if (p.phi < 0) p.phi += 2 * M_PI;
    p.phi -= M_PI;
}