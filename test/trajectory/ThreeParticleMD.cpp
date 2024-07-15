#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/Particle.h"
#include "lettuce/Utilities.h"
#include "lettuce/Integration.h"
#include "lettuce/LennardJones.h"
#include "lettuce/CircleDistribution.h"
#include "lettuce/CirclesIntersectionFuncs.h"


int main(){
    const double sigma = 1;
    const double epsilon = 1;
    const double mass = 1;
    const double radius = 1;
    const double cutoff = 20;
    const double dt = 0.01;
    const double Time = 1000.;
    const double stepT = 10;

    Species<double> species1 {mass, radius};
    Species<double> species2 {1.2 * mass, 0.9 * radius};

    std::vector <Species<double>> allSpecies {species1, species2};

    Particle<double> P1{{{0.0, 0.0}}, {{0.0, 0.0}}, 0};
    Particle<double> P2{{{1.2, 0.0}}, {{0.0, 0.0}}, 0};
    Particle<double> P3{{{0.0, 1.2}}, {{0.0, 0.0}}, 1};
    Particle<double> P4{{{1.2, 1.2}}, {{0.0, 0.0}}, 0};

    std::vector<Particle<double>> particles {P1, P2, P3, P4};

    LennardJonesForce<double> LJForce(epsilon, sigma, cutoff);
    LennardJonesPotential<double> LJPotential(epsilon, sigma, cutoff);

    VelocityVerletIntegrator<double, LennardJonesForce<double>> Method;

    std::ofstream positionFile ("PosMD.dat");
    std::ofstream kineticEnergyFile ("KineticEnergyMD.dat");
    std::ofstream PotentialEnergyFile ("PotentialEnergyMD.dat");

    int numSteps = static_cast<int>(Time / stepT);

    for (int i = 0; i < numSteps; ++i){
        integrate(particles, allSpecies, dt, Time, std::move(LJForce), Method);
        writePositionToFile(particles, positionFile);
        writeKineticEToFile(particles, allSpecies, kineticEnergyFile);
        writePotentialEToFile(particles, allSpecies, LJPotential, PotentialEnergyFile);     //why it does not need std::move()?
    }
    
    return 0;
}