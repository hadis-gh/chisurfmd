#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include "Particle.h"

enum class InitialParticlesConfiguration {
    RANDOM_CIRCLES,
    DLA
};

template<typename T>
Particle<T> startParticleRandom (const T& radius, const T& areaWidth, std::mt19937 &gen){   
    Particle<T> newParticle;

    std::uniform_real_distribution<> randomPos(0, areaWidth);
    std::uniform_int_distribution<> randomSide(1,4); 

    if (randomSide.operator()(gen) == 1) {
        newParticle.r[0] = randomPos(gen);
        newParticle.r[1] = 0;
    }
    else if(randomSide.operator()(gen) == 2){
        newParticle.r[0] = 0;
        newParticle.r[1] = randomPos(gen);
    }
    else if(randomSide.operator()(gen) == 3){
        newParticle.r[0] = areaWidth;
        newParticle.r[1] = randomPos(gen);        
    }
    else {
        newParticle.r[0] = randomPos(gen);
        newParticle.r[1] = areaWidth;
    }

    return newParticle;
}

template<typename T>
Vec<T> shootToCenter(const Particle<T> &startParticle, const T &areaWidth) {
    return -startParticle.r + areaWidth / static_cast<T>(2);
}

template <typename T>
bool hasOverlap (const Particle<T> &p1,const Particle<T> &p2, const T &radius){
    T distance2 = (p1.r - p2.r).abs2();
    return distance2 < (2 * radius);
}

template<typename T>
Particle<T> placeRandomParticle(const std::vector <Particle<T>> &particles, std::uniform_real_distribution<> dis, const T &radius, std::mt19937 &gen)
{
    Particle<T> newParticle;
    newParticle.r[0] = dis(gen);
    newParticle.r[1] = dis(gen);

    bool overlap = false;
    
    for (auto p: particles){
        if(hasOverlap(newParticle, p, radius)){
            overlap = true;
            break;
        }
    }
    if(overlap)      
        newParticle.r = NAN;
    return newParticle;
}

template <typename T>
std::vector <Particle<T>> distRandomParticles (const int &particlesNum, const T &radius, const T &L, std::mt19937 &gen){
    std::vector <Particle<T>> particles;
    std::uniform_real_distribution<> dis(0, L);
    
    while (particles.size() < particlesNum){
        auto newParticle = placeRandomParticle(particles, dis, radius, gen);
        if(!std::isnan(newParticle.r[0]))
           particles.push_back(newParticle);
    }
    return particles;
}

template<typename T>
std::vector<Particle<T>> distParticlesDLA (const int &shootNum, const T &radius, const T &L, std::mt19937 &gen){
    Particle<T> target = {{{L / 2, L / 2}}};
    std::vector<Particle<T>> finalParticles;
    finalParticles.push_back(target);

    for (int j = 0; j < shootNum; ++j) {
        Particle<T> newParticle = startParticleRandom(radius, L, gen);
        Vec<T> direction = shootToCenter(newParticle, L);

        Particle<T> endPoint = findStopPointAll(newParticle, radius, direction, finalParticles);
        if (!std::isnan(endPoint.r[0])){
            finalParticles.push_back(endPoint);
        }
    }
    return finalParticles;
}

template<typename T>
std::vector<Particle<T>> initialParticles (const int &particlesNum, const std::vector<Species<T>> &allSpecies, int &speciesNum, const T &L, std::mt19937 &gen, const std::string& configuation){
    std::vector<Particle<T>> particles;
    if (configuation == "RANDOM"){
        particles = distRandomParticles (particlesNum, allSpecies[speciesNum].radius, L, gen);
    }else if(configuation == "DLA"){
        particles = distParticlesDLA (particlesNum, allSpecies[speciesNum].radius, L, gen);
    }else if(configuation == "THREE"){
        particles.reserve(4);
        particles.push_back({{{0.0, 0.0}}, {{0.0, 0.0}}, 0});
        particles.push_back({{{1.2, 0.0}}, {{0.0, 0.0}}, 0});
        particles.push_back({{{0.0, 1.2}}, {{0.0, 0.0}}, 1});
        particles.push_back({{{1.2, 1.2}}, {{0.0, 0.0}}, 0});

        return particles;
    }
    for (auto &p: particles){
        p.species = speciesNum;
    }
    return particles;
}
