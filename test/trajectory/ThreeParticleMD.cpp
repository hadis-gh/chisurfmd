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

const auto sigma6 = pow(sigma,6);
const auto sigma12 = sigma6 * sigma6;

template<typename T>
struct Particle{
    Vec<T> r;
    Vec<T> v;
    T mass;
};

template<typename T>
Vec<T> calForceTwo(const Particle<T> &p1, const Particle<T> &p2){
    T r = (p2.r - p1.r).abs2();
    if (r == 0) return {{0, 0}};
    auto rPower6 = pow(r, 6);
    
    T ljForce = 48.0 * epsilon * (sigma6 * sigma6 /(rPower6 * rPower6 * r) - 0.5 * sigma6 /(rPower6 * r));
    return ljForce * (p2.r - p1.r)/r;
}

template<typename T>
T ljPotential(T &r){
    if (r == 0) return 0;
    T r6 = pow(r, 6);
    T r12 = r6 * r6;
    return 4.0 * epsilon * (sigma12 /(r12) - 0.5 * sigma6 /(r6));
}

template<typename T>
T ljForce(T &r){
    if (r == 0) return 0;
    T r6 = pow(r, 6);
    T r12 = r6 * r6;
    return 48.0 * epsilon * (sigma12 /(r12 * r) - 0.5 * sigma6 /(r6 * r));
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
std::vector<Vec<T>> calculateAccelerations(std::vector<Particle<T>> &particles) {
    std::vector<Vec<T>> accelarations;
    for (auto &p : particles) {
        accelarations.push_back(calAccelaration(p, particles) / p.mass);
    }
    return accelarations;
}

template<typename T>        //Euler Integration
void updateStateEuler(std::vector<Particle<T>> &particles, const T &dt){
    for (auto &p : particles){
        p.v += calAccelaration(p, particles) * dt;
        p.r += p.v * dt;
    }
}

template<typename T>        //Velocity Vernel Integration
void updateStateVV(std::vector<Particle<T>> &particles, const T &dt) {
    for (auto &p : particles) {
        p.r += p.v * dt + calAccelaration(p, particles)/2 * dt * dt;
    }
    std::vector<Vec<T>> old_accelerations;
    for (const auto &p : particles) {
        old_accelerations.push_back(calAccelaration(p, particles));
    }
    std::vector<Vec<T>> accelarations = calculateAccelerations(particles);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].v += (old_accelerations[i] + accelarations[i])/2 * dt;
    }
}

template<typename T>
void writeParticleProperties(const std::vector<Particle<T>> &particles) {
    for (const auto &p : particles) {
        std::cout << "Position: " << p.r << "\t"
                  << "Velocity: " << p.v << "\t" << std::endl;
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

template<typename T>
std::vector<T> linspace(const T& start, const T& end, const int& points) {
    std::vector<T> result;
    T step = (end - start) / (points - 1);
    for (int i = 0; i < points; ++i) {
        result.push_back(start + i * step);
    }
    return result;
}

template<typename T>
void writePlotData(const T &startR, const T &endR, const int &numSpace){
    std::vector<double> distances = linspace(startR, endR, numSpace);
    std::vector<double> ljForceValue;
    std::vector<double> ljPotentialValue;

    for (auto &r: distances){
        ljForceValue.push_back(ljForce(r));
        ljPotentialValue.push_back(ljPotential(r));
    }

    std::ofstream ljForceFile ("plotLJforce.txt");
    for (int i= 0; i < distances.size(); ++i){
        ljForceFile << distances[i] << " " << ljForceValue[i] << std::endl;
    }

    std::ofstream ljPotentialFile ("plotLJpotential.txt");
    for (int i= 0; i < distances.size(); ++i){
        ljPotentialFile << distances[i] << " " << ljPotentialValue[i] << std::endl;
    }
}

int main(){
    Particle<double> P1{{{0.0, 0.0}}, {{0.0, 0.0}}, mass};
    Particle<double> P2{{{10e-5, 0.0}}, {{0.0, 0.0}}, mass};
    Particle<double> P3{{{0.0, 10e-5}}, {{0.0, 0.0}}, mass};

    std::vector<Particle<double>> particles {P1, P2, P3};

    int numSteps = 10;

    std::ofstream file ("particlesPosMD.txt");

    for (int i = 0; i < numSteps; ++i){
        updateStateEuler(particles, dt);
        writeToFile(particles, file);
    }

    const double startR = 0.9 * sigma;
    const double endR = 3.0 * sigma;
    const int numSpace = 400;

    writePlotData(startR, endR, numSpace);

    return 0;
}