#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"

template<typename T>
struct Species
{
    T mass;
    T momentOfInertia;
    T radius;

    Species(T mass, T momentOfInertia, T radius)
        : mass(mass), momentOfInertia(momentOfInertia), radius(radius)
    {}
};

template<typename T>
struct ParticleDot
{
    using value_type = T;

    Vec<value_type> r, v;
    int h,d;
    
    unsigned int species;

    ParticleDot() : species(0), r({0, 0}), v({0, 0}), h(1), d(1) {}

    ParticleDot(unsigned int speciesIndex) : species(speciesIndex), r({0, 0}), v({0, 0}) {}

    ParticleDot(unsigned int speciesIndex, Vec<value_type> position) : species(speciesIndex), r(position), v({0, 0}) {}

    ParticleDot(unsigned int speciesIndex, Vec<value_type> position, Vec<value_type> velocity)
        : species(speciesIndex), r(position), v(velocity) {}
};

template<typename Particle, typename SFINAE=void>
struct CreateRandomParticle;

template<typename Particle>
auto createRandomParticle(std::mt19937& gen){
    return CreateRandomParticle<Particle>::createRandomParticle(gen);
}

template<typename T>
struct CreateRandomParticle<ParticleDot<T>>
{
    static ParticleDot<T> createRandomParticle(std::mt19937& gen)
    {
        std::uniform_real_distribution<T> randomPos(0, 1);
        ParticleDot<T> p;
        p.r[0] = randomPos(gen);
        p.r[1] = randomPos(gen);
        return p;
    }
};

template<typename Particle, typename SFINAE=void>
struct CreateTwoParticle;

template<typename Particle, typename T>
auto createTwoParticle(std::mt19937& gen, const T& areaL) {
    return CreateTwoParticle<Particle>::createTwoParticle(gen, areaL);
}

template<typename T>
struct CreateTwoParticle<ParticleDot<T>> {
    static std::vector<ParticleDot<T>> createTwoParticle(std::mt19937& gen, const T& areaL) { // pass potential as a parameter with equilibr distance
        ParticleDot<T> p1, p2;
        
        T eqDis = 1.123;
        p1.r[0] = areaL / 2;
        p1.r[1] = areaL / 2;
        
        p2.r[0] = areaL / 2 + eqDis;
        p2.r[1] = areaL / 2 + eqDis;
        
        return std::vector<ParticleDot<T>>{p1, p2};
    }
};

template<typename Particle, typename SFINAE=void>
struct DegreesOfFreedom;

template<typename Particle>
constexpr auto degreesOfFreedom(){
    return DegreesOfFreedom<std::decay_t<Particle>>::degreesOfFreedom();
}

template<typename T>
struct DegreesOfFreedom<ParticleDot<T>>
{
    static constexpr int degreesOfFreedom()
    {
        return 2;
    }
};

// getting generalized positions
template<typename Particle, typename SFINAE=void>
struct GetGeneralizedPositions;

template<typename Particle>
auto getGeneralizedPositions(Particle&& p){
    return GetGeneralizedPositions<std::decay_t<Particle>>::getGeneralizedPositions(std::forward<Particle>(p));
}

template<typename T>
struct GetGeneralizedPositions<ParticleDot<T>>
{
    static const Vec<T,2>& getGeneralizedPositions(const ParticleDot<T>& p)
    {
        return p.r;
    }
};

// setting generalized positions
template<typename Particle, typename SFINAE = void>
struct SetGeneralizedPositions;

template<typename Particle>
auto setGeneralizedPositions(Particle&& p, const auto& new_positions) {
    return SetGeneralizedPositions<std::decay_t<Particle>>::setGeneralizedPositions(std::forward<Particle>(p), new_positions);
}

template<typename T>
struct SetGeneralizedPositions<ParticleDot<T>>
{
    static void setGeneralizedPositions(ParticleDot<T>& p, const Vec<T, 2>& new_positions)
    {
        p.r = new_positions;
    }
};

// getting generalized velocities
template<typename Particle, typename SFINAE = void>
struct GetGeneralizedVelocities;

template<typename Particle>
auto getGeneralizedVelocities(Particle&& p)
{
    return GetGeneralizedVelocities<std::decay_t<Particle>>::getGeneralizedVelocities(std::forward<Particle>(p));
}

template<typename T>
struct GetGeneralizedVelocities<ParticleDot<T>>
{
    static const Vec<T, 2>& getGeneralizedVelocities(const ParticleDot<T>& p)
    {
        return p.v;
    }
};

// setting generalized velocities
template<typename Particle, typename SFINAE = void>
struct SetGeneralizedVelocities;

template<typename Particle>
auto setGeneralizedVelocities(Particle&& p, const auto& new_velocities)
{
    return SetGeneralizedVelocities<std::decay_t<Particle>>::setGeneralizedVelocities(std::forward<Particle>(p), new_velocities);
}

template<typename T>
struct SetGeneralizedVelocities<ParticleDot<T>>
{
    static void setGeneralizedVelocities(ParticleDot<T>& p, const Vec<T, 2>& new_velocities)
    {
        p.v = new_velocities;
    }
};

//Force calculations
template<typename T, typename Force>
Vec<T> calForceTwo(const ParticleDot<T> &p1, const ParticleDot<T> &p2, const T& boxPBC, Force &&force){
    Vec<T> dr = p2.r - p1.r;
    
    for (int i = 0; i < 2; ++i) {
        if (dr[i] > boxPBC / 2) { dr[i] -= boxPBC; }
        else if (dr[i] < -boxPBC / 2) { dr[i] += boxPBC; }
    }
    
    T r = dr.abs();
    if (r == 0) return {{0, 0}};
    
    T f = force(p1, p2, dr, r);
    return f * dr / r;
}

template<typename T, typename Force>
Vec<T> calTotalForce(const ParticleDot<T> &p1, const std::vector<ParticleDot<T>> &particles, const T& boxPBC, Force &&force){
    Vec<T> f;
    for (auto &p : particles){
        if (p1.r != p.r){
            f += calForceTwo(p, p1, boxPBC, std::forward<Force>(force));
        }
    }
    return f;
}

template<typename T, typename Force>
Vec<T> calAccelaration(const ParticleDot<T> &p1, const std::vector<ParticleDot<T>> &particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Force &&force){
    T mass = allSpecies[p1.species].mass;
    return calTotalForce(p1, particles, boxPBC, std::forward<Force>(force))/ mass;
}

template<typename TParticle, typename Force, typename T = typename TParticle::value_type>
std::vector<Vec<T>> calAllAccelerations(const std::vector<TParticle> &particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Force &&force) {
    std::vector<Vec<T, degreesOfFreedom<TParticle>()>> accelerations(particles.size());
    for (unsigned int i = 0; i < particles.size(); ++i) {
        accelerations[i] = calAccelaration(particles[i], particles, allSpecies, boxPBC, std::forward<Force>(force));
    }
    return accelerations;
}

template<typename T>
void implementPBC(ParticleDot<T>& p, const T& boxPBC) {    
    for (int i = 0; i < 2; ++i) {
        if (p.r[i] > boxPBC) { 
            p.r[i] -= boxPBC; 
        } else if (p.r[i] < 0) { 
            p.r[i] += boxPBC; 
        }
    }
}