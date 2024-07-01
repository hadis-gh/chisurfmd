#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/CircleDistribution.h"
#include "lettuce/CirclesIntersectionFuncs.h"

// Derivative is based on Lennar Jones Potential: 4.0 * epsilon * (pow(sigma / r, 12) - pow(sigma / r, 6))

const float sigma = 1e-10;
const float epsilon = 1e-21;
const float mass = 1e-25;

template<typename T>
struct Particle{
    Vec<T> r;
    Vec<T> v;
    Vec<T> a;
    T mass;
};

template<typename T>
Vec<T> calForceTwo(const Particle<T> &p1, const Particle<T> &p2){
    T r = (p2.r - p1.r).abs2();
    return -48.0 * epsilon * (pow(sigma / r, 13) - 0.5 * pow(sigma / r, 7));
}

template<typename T>
Vec<T> calTotalForce(const Particle<T> &p1, const std::vector<Particle<T>> &particles){
    Vec<T> force;
    for (auto &p : particles){
        if (p1.r != p.r){
            force += calForceTwo(p1, p);
        }
    }
    return force;
}

template<typename T>
Vec<T> calAccelaration(const Particle<T> &p1, const std::vector<Particle<T>> &particles){
    return calTotalForce(p1, particles)/ mass;
}

template<typename T>
void updateState(Particle<T> &p1, const std::vector<Particle<T>> &particles, const T &dt){
    p1.a = calAccelaration(p1, particles);
    p1.v += p1.a * dt;
    p1.r += p1.v * dt;
}

template<typename T>
void writeParticleProperties(const Particle<T> &p){
    std::cout << "Position: \t"       << p.r << std::endl;
    std::cout << "Velocity: \t"       << p.v << std::endl;
    std::cout << "Acceleration: \t"   << p.a << std::endl;
}

int main(){
    Particle<float> P1{{{0.0, 0.0}}, {{1e-5, 1e-5}}};
    Particle<float> P2{{{5e-10, 0.0}}, {{0.0, 0.0}}};
    Particle<float> P3{{{0.0, 5e-10}}, {{0.0, 0.0}}};

    const float dt = 1e-20;

    std::vector<Particle<float>> particles {P1, P2, P3};

    std::cout << "Force of P2: \t"      << calForceTwo(P1, P2) << std::endl;

    std::cout <<"Before Update:-----------------" << std::endl;
    writeParticleProperties(P1);
    std::cout <<"\nAfter Update:-----------------" << std::endl;
    updateState(P1, particles, dt);
    writeParticleProperties(P1);
    return 0;
}
