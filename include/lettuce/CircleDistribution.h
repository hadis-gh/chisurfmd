#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include "Circle.h"
#include "lettuce/Particle.h"
#include "lettuce/CirclesIntersectionFuncs.h"

template<typename T>
Circle<T> startCircleRandom (const T& radius, const T& areaWidth, std::mt19937 &gen){   
    Circle<T> newCircle;
    newCircle.r = radius;

    std::uniform_real_distribution<> randomPos(0, areaWidth);
    std::uniform_int_distribution<> randomSide(1,4); 

    if (randomSide.operator()(gen) == 1) {
        newCircle.c[0] = randomPos(gen);        //why can't we write newCircle.c={{randomPos(gen), 0}};
        newCircle.c[1] = 0;
    }
    else if(randomSide.operator()(gen) == 2){
        newCircle.c[0] = 0;
        newCircle.c[1] = randomPos(gen);
    }
    else if(randomSide.operator()(gen) == 3){
        newCircle.c[0] = areaWidth;
        newCircle.c[1] = randomPos(gen);        
    }
    else {
        newCircle.c[0] = randomPos(gen);
        newCircle.c[1] = areaWidth;
    }

    return newCircle;
}

template<typename T>
Vec<T> shootToCenter(const Circle<T> &startCircle, const T &areaWidth) {
    return -startCircle.c + areaWidth / static_cast<T>(2);
}

template <typename T>
bool has_overlap (const Circle<T> &c1,const Circle<T> &c2){
    T distance2 = (c1.c - c2.c).abs2();
    const auto radiuses (c1.r+c2.r);
    return distance2 < radiuses*radiuses;
}

template<typename T>
Circle<T> placeRandomCircle(const std::vector <Circle<T>> &circles, std::uniform_real_distribution<> dis, const T &radius, std::mt19937 &gen)
{
    Circle<T> newCircle;
    newCircle.c[0] = dis(gen);
    newCircle.c[1] = dis(gen);
    newCircle.r = radius;

    bool overlap = false;
    
    for (auto circle: circles){
        if(has_overlap(newCircle, circle)){
            overlap = true;
            break;
        }
    }
    if(overlap)      
        newCircle.c = NAN;
    return newCircle;
}

template <typename T>
std::vector <Circle<T>> distRandomCircles (const int &circlesNum, const T &L, const T &radius, std::mt19937 &gen){
    std::vector <Circle<T>> circles;
    std::uniform_real_distribution<> dis(0, L);
    
    while (circles.size() < circlesNum){
        auto newCircle = placeRandomCircle(circles, dis, radius, gen);
        if(!std::isnan(newCircle.c[0]))
           circles.push_back(newCircle);
    }
    return circles;
}

template<typename T>
std::vector<Particle<T>> distRandomParticles (const int &particlesNum, const T &L, const T &radius, std::mt19937 &gen){
    std::vector<Particle<T>> randomParticles(particlesNum);
    std::vector<Circle<T>> randomCircles(particlesNum);
    randomCircles = distRandomCircles(particlesNum, L, radius, gen);
    for (int i = 0; i < particlesNum; ++i){
        randomParticles[i].r = randomCircles[i].c;
    }
    return randomParticles;
}

template<typename T>
std::vector<Circle<T>> distCirclesDLA (const int &shootNum, const T &L, const T &radius, std::mt19937 &gen){
    Circle<T> target = {{{L / 2, L / 2}}, radius};
    std::vector<Circle<T>> finalCircles;
    finalCircles.push_back(target);

    for (int j = 0; j < shootNum; ++j) {
        Circle<T> newCircle = startCircleRandom(radius, L, gen);
        Vec<T> direction = shootToCenter(newCircle, L);

        Circle<T> endPoint = findStopPointAll(newCircle, direction, finalCircles);
        if (!std::isnan(endPoint.c[0])){
            finalCircles.push_back(endPoint);
        }
    }
    return finalCircles;
}

template<typename T>
std::vector<Particle<T>> distParticleDLA (const int &shootNum, const T &L, const T &radius, std::mt19937 &gen){
    std::vector<Circle<T>> finalCircles = distCirclesDLA (shootNum, L, radius, gen);
    std::vector<Particle<T>> finalParticles(finalCircles.size());

    for (int i =0; i < finalCircles.size(); ++i){
        finalParticles[i].r = finalCircles[i].c;
    }
    return finalParticles;
}

template<typename T>
void writeCircles(T begin, T end, const std::string& fname){
        std::ofstream output_file(fname);
        for (auto c = begin; c != end; ++c)
            output_file << c->c << ", " << c->r <<"\n";
}

template<typename T>
void writeParticle(const std::vector<Particle<T>> &particles, T radius, const std::string& fname){
        std::ofstream output_file(fname);
        for (auto p: particles)
            output_file << p.r[0] << ", " << p.r[1] << ", " << radius <<"\n";
}

template<typename T>
std::vector<Particle<T>> initialParticles(const unsigned int &particlesNum, const std::vector<Species<T>> &allSpecies, int &speciesNum, const T &L, std::mt19937 &gen, const std::string& configuration) {
    std::vector<Particle<T>> particles;

    std::cout << "Configuration: " << configuration << std::endl;

    if (configuration == "RANDOM") {
        particles = distRandomParticles(particlesNum, allSpecies[speciesNum].radius, L, gen);
    } else if (configuration == "DLA") {
        particles = distParticleDLA(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } else if (configuration == "THREE") {
        particles.reserve(4);
        particles.push_back({{{0.0, 0.0}}, {{0.0, 0.0}}, 0});
        particles.push_back({{{1.2, 0.0}}, {{0.0, 0.0}}, 0});
        particles.push_back({{{0.0, 1.5}}, {{0.0, 0.0}}, 0});
        particles.push_back({{{1.2, 1.2}}, {{0.0, 0.0}}, 0});
    } else if (configuration == "FILE"){

        std::ifstream configutation("/home/hadis/custom_vector/build/test/configuration.dat");
        int numParticles = 0;
        std::string line;
        while (std::getline(configutation, line)) {
            ++numParticles;
        }
        particles.reserve(numParticles);
        for (auto& p : particles) {
            configutation >> p.r[0] >> p.r[1] >> p.v[0] >> p.v[1];
        }
    }
    for (auto &p: particles){
        p.species = speciesNum;
    }
    return particles;
}

