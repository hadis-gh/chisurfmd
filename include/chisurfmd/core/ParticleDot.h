#pragma once

#include <random>
#include <vector>
#include <cmath>
#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/Species.h"

template<typename T>
struct ParticleDot
{
    using value_type = T;

    // Todo r,v,h,d is not good naming.
    Vec<value_type> position, velocity;
    int8_t handedness, alignment;
    bool fixed;

    unsigned int species;

    ParticleDot()
        : species(0), position({0, 0}), velocity({0, 0}), handedness(1), alignment(1), fixed(false) {}

    ParticleDot(unsigned int speciesIndex)
        : species(speciesIndex), position({0, 0}), velocity({0, 0}), handedness(1), alignment(1), fixed(false) {}

    ParticleDot(unsigned int speciesIndex, Vec<value_type> position)
        : species(speciesIndex), position(position), velocity({0, 0}), handedness(1), alignment(1), fixed(false) {}

    ParticleDot(unsigned int speciesIndex, Vec<value_type> position, Vec<value_type> velocity)
        : species(speciesIndex), position(position), velocity(velocity), handedness(1), alignment(1), fixed(false) {}
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
        p.position[0] = randomPos(gen);
        p.position[1] = randomPos(gen);
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
        p1.position[0] = areaL / 2;
        p1.position[1] = areaL / 2;
        
        p2.position[0] = areaL / 2 + eqDis;
        p2.position[1] = areaL / 2 + eqDis;
        
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
        return p.position;
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
        p.position = new_positions;
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
        return p.velocity;
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
        p.velocity = new_velocities;
    }
};

//Force calculations
template<typename T, typename Force>
Vec<T, 2> calForceTwo(const ParticleDot<T>& p1, const ParticleDot<T>& p2, const T& boxPBC, Force&& force) {
    Vec<T, 2> dr = p2.position - p1.position;
    
    for (int i = 0; i < 2; ++i) {
        if (dr[i] > boxPBC / 2) dr[i] -= boxPBC;
        else if (dr[i] < -boxPBC / 2) dr[i] += boxPBC;
    }
    
    const T r = dr.abs();
    if (r <= std::numeric_limits<T>::epsilon()) return {{0, 0}};
    
    return force(p1, p2, dr, r);
}

template<typename T, typename Force>
Vec<T, 2> calTotalForce(const ParticleDot<T>& p1, const std::vector<ParticleDot<T>>& particles, 
                         const T& boxPBC, Force&& force) {
    Vec<T, 2> totalForce{{0, 0}};
    for (const auto& p : particles) {
        if (&p1 != &p) {
            totalForce += calForceTwo(p, p1, boxPBC, std::forward<Force>(force));
        }
    }
    return totalForce;
}

template<typename T, typename Force>
Vec<T, 2> calAcceleration(const ParticleDot<T>& p1, const std::vector<ParticleDot<T>>& particles,
                           const std::vector<Species<T>>& allSpecies, const T& boxPBC, Force&& force) {
    const T mass = allSpecies[p1.species].mass;
    return calTotalForce(p1, particles, boxPBC, std::forward<Force>(force)) / mass;
}

template<typename TParticle, typename Force, typename T = typename TParticle::value_type>
std::vector<Vec<T, 2>> calAllAccelerations(const std::vector<TParticle>& particles, 
                                           const std::vector<Species<T>>& allSpecies, 
                                           const T& boxPBC, Force&& force) {
    std::vector<Vec<T, 2>> accelerations;
    accelerations.reserve(particles.size());
    
    for (const auto& p : particles) {
        accelerations.push_back(
            calAcceleration(p, particles, allSpecies, boxPBC, std::forward<Force>(force)));
    }
    return accelerations;
}

template<typename T>
void implementPBC(ParticleDot<T>& p, const T& boxPBC) {    
    for (int i = 0; i < 2; ++i) {
        p.position[i] = std::fmod(p.position[i], boxPBC);
        if (p.position[i] < 0) p.position[i] += boxPBC;
    }
}