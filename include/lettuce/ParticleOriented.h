#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/ParticleDot.h"

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
Vec<T, 3> calForceTwo(const ParticleOriented<T>& p1, const ParticleOriented<T>& p2, const T& boxPBC, Force&& force) {
    Vec<T, 2> dr = p2.r - p1.r;
    T deltaPhi = p2.phi - p1.phi;

    for (int i = 0; i < 2; ++i) {   //make extra func
        if (dr[i] > boxPBC / 2) { dr[i] -= boxPBC; }
        else if (dr[i] < -boxPBC / 2) { dr[i] += boxPBC; }
    }
    
    T r = dr.abs();
    if (r == 0) return {{0, 0}};
    
    Vec<T, 2> f = force(p1, p2, dr, r);
    Vec<T, 3> force_vec;

    force_vec[0] = f[0] * dr[0] / r;
    force_vec[1] = f[0] * dr[1] / r;
    force_vec[2] = f[1];
    
    return force_vec;
}

template<typename T, typename Force>
Vec<T, 3> calTotalForce(const ParticleOriented<T>& p1, const std::vector<ParticleOriented<T>>& particles, const T& boxPBC, Force&& force) {
    Vec<T, 3> total_force;
    for (const auto& p : particles) {
        if (p1.r != p.r) {
            total_force += calForceTwo(p, p1, boxPBC, std::forward<Force>(force));
        }
    }
    return total_force;
}

template<typename T, typename Force>
Vec<T, 3> calAcceleration(const ParticleOriented<T>& p1, const std::vector<ParticleOriented<T>>& particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Force&& force) {
    T mass = allSpecies[p1.species].mass;
    T MI = allSpecies[p1.species].momentOfInertia;

    Vec<T, 3> total_force = calTotalForce(p1, particles, boxPBC, std::forward<Force>(force));
    Vec<T, 2> acceleration = {{total_force[0] / mass, total_force[1] / mass}};
    T angular_acceleration = total_force[2] / MI;

    return {{acceleration[0], acceleration[1], angular_acceleration}};
}

template<typename T, typename Force>
std::vector<Vec<T, 3>> calAllAccelerations(const std::vector<ParticleOriented<T>>& particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Force&& force) {
    std::vector<Vec<T, 3>> accelerations(particles.size());
    for (unsigned int i = 0; i < particles.size(); ++i) {
        accelerations[i] = calAcceleration(particles[i], particles, allSpecies, boxPBC, std::forward<Force>(force));
    }
    return accelerations;
}

template<typename T>
void implementPBC(ParticleOriented<T>& p, const T& boxPBC) {
    implementPBC(static_cast<ParticleDot<T>&>(p), boxPBC);

    T phi = p.phi;
    phi = std::fmod(phi + M_PI, 2 * M_PI); 
    if (phi < 0) phi += 2 * M_PI;
    p.phi = phi - M_PI;
}