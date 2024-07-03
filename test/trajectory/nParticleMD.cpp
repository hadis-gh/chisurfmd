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

const double sigma = 1e-10;
const double epsilon = 1e-21;
const double mass = 1e-25;
const double dt = 1e-20;

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
    if (r == 0) return {{0, 0}};
    T ljForce = 48.0 * epsilon * (pow(sigma, 12)/pow(r, 13) - 0.5 * pow(sigma,6)/pow(r, 7));
    return ljForce * (p2.r - p1.r)/r;
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
    return calTotalForce(p1, particles)/ p1.mass;
}
template<typename T>
void calculateAccelerations(std::vector<Particle<T>> &particles) {
    for (auto &p : particles) {
        p.a = calTotalForce(p, particles) / p.mass;
    }
}

template<typename T>        //Euler Integration
void updateStateEuler(std::vector<Particle<T>> &particles, const T &dt){
    for (auto &p : particles){
        p.v += p.a * dt;
        p.r += p.v * dt;
        p.a = calAccelaration(p, particles);
    }
}

template<typename T>        //Velocity Vernel Integration
void updateStateVV(std::vector<Particle<T>> &particles, const T &dt) {
    for (auto &p : particles) {
        p.r += p.v * dt + p.a/2 * dt * dt;
    }
    std::vector<Vec<T>> old_accelerations;
    for (const auto &p : particles) {
        old_accelerations.push_back(p.a);
    }
    calculateAccelerations(particles);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].v += (old_accelerations[i] + particles[i].a)/2 * dt;
    }
}

template<typename T>
void writeParticleProperties(const std::vector<Particle<T>> &particles) {
    for (const auto &p : particles) {
        std::cout << "Position: " << p.r << "\t"
                  << "Velocity: " << p.v << "\t"
                  << "Acceleration: " << p.a << std::endl;
    }
}

template<typename T>
void debuggingPrintInfo(std::vector<Particle<T>> &particles, const int &numSteps, const T &dt){
    std::cout <<"Initial states:" << std::endl;
    writeParticleProperties(particles);

    for (int step = 0; step < numSteps; ++step) {
        updateStateVV(particles, dt);

        std::cout << "\nAfter Update step " << step + 1 << ":" << std::endl;
        writeParticleProperties(particles);
    }
}

template<typename T>
void writeToFile(const std::vector<Particle<T>> &particles, std::ofstream &file){
    for (auto &p: particles){
        file << p.r[0] << " " << p.r[1] << " ";
    }
    file << "\n";
}

int main(){
    Particle<double> P1{{{0.0, 0.0}}, {{0.0, 0.0}}, {{0.0, 0.0}}, mass};
    Particle<double> P2{{{10e-9, 0.0}}, {{0.0, 0.0}}, {{0.0, 0.0}}, mass};
    Particle<double> P3{{{0.0, 10e-9}}, {{0.0, 0.0}}, {{0.0, 0.0}}, mass};

    std::vector<Particle<double>> particles {P1, P2, P3};
    int numSteps = 10;

    std::ofstream file ("particlesPosMD.txt");

    for (int i = 0; i < numSteps; ++i){
        updateStateEuler(particles, dt);
        writeToFile(particles, file);
    }

    return 0;
}