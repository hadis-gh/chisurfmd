#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include "Circle.h"
#include "lettuce/Particle.h"


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