#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include "Circle.h"
#include "lettuce/Particle.h"

enum InitialParticlesConfiguration {
    RANDOM_CIRCLES,
    DLA
};

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
std::vector<Circle<T>> distCirclesPBC (const int &circlesNum, const T &L, const T &radius, std::mt19937 &gen){
    std::vector<Circle<T>> repCircles = distRandomCircles (circlesNum, L, radius, gen);
    
    for (const auto &circle : repCircles) {
        Vec<T> c = circle.c;
        T r = circle.r;
        // left, right, bottom, up
        if (c[0] - r < 0)
            repCircles.emplace_back(Vec<T>({c[0] + L, c[1]}), r);
        if (c[0] + r >= L) 
            repCircles.emplace_back(Vec<T>({c[0] - L, c[1]}), r);
        if (c[1] - r < 0) 
            repCircles.emplace_back(Vec<T>({c[0], c[1] + L}), r);
        if (c[1] + r >= L)
            repCircles.emplace_back(Vec<T>({c[0], c[1] - L}), r);
        //corners
        if (c[0] - r < 0 && c[1] - r < 0)
            repCircles.emplace_back(Vec<T>({c[0] + L, c[1] + L}), r);
        if (c[0] + r >= L && c[1] - r < 0) 
            repCircles.emplace_back(Vec<T>({c[0] - L, c[1] + L}), r);
        if (c[0] - r < 0 && c[1] + r >= L)
            repCircles.emplace_back(Vec<T>({c[0] + L, c[1] - L}), r);
        if (c[0] + r >= L && c[1] + r >= L)
            repCircles.emplace_back(Vec<T>({c[0] - L, c[1] - L}), r);
    }
    return repCircles;
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
    std::vector<Circle<T>> finalCircles;
    std::vector<Particle<T>> finalParticles;

    finalCircles = distCirclesDLA (shootNum, L, radius, gen);

    for (int i =0; i < finalCircles.size(); ++i){
        finalParticles[i].r = finalCircles[i].c;
    }                                               //can we write just because of Vec<T>: `finalParticles.r = finalCircles.c;`

    return finalParticles;
}

template<typename T>
void writeCircles(T begin, T end, const std::string& fname){
        std::ofstream output_file(fname);
        for (auto c = begin; c != end; ++c)
            output_file << c->c << ", " << c->r <<"\n";
}


template<typename T>
std::vector<Particle<T>> initialParticles (const int &particlesNum, const T &L, const T &radius, std::mt19937 &gen, InitialParticlesConfiguration configuation){
    std::vector<Particle<T>> particles(particlesNum);

    if (configuation == RANDOM_CIRCLES){
        particles = distRandomParticles (particlesNum, L, radius, gen);
    }else if(configuation == DLA){
        particles = distParticleDLA (particlesNum, L, radius, gen);
    }
    return particles;
}