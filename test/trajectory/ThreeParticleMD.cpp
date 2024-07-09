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
    const double cutoff = 20;
    const double dt = 0.01;
    const double Time = 1000.;

    IntegrationMethod integrationMethod = EULER; 

    Particle<double> P1{{{0.0, 0.0}}, {{0.0, 0.0}}, mass};
    Particle<double> P2{{{1.2, 0.0}}, {{0.0, 0.0}}, mass};
    Particle<double> P3{{{0.0, 1.2}}, {{0.0, 0.0}}, mass};
    Particle<double> P4{{{0.0, 5.1}}, {{0.0, 0.0}}, mass};

    std::vector<Particle<double>> particles {P1, P2, P3, P4};

    LennardJonesForce<double> LJForce(epsilon, sigma, cutoff);
    LennardJonesPotential<double> LJPotential(epsilon, sigma, cutoff);

    integrate(particles, dt, Time, LJForce, integrationMethod);

    std::vector<double> kineticEnergy;
    kineticEnergy = calculateKineticEnergy(particles, Time, dt);

    std::vector<double> potentialEnergyLJ;
    potentialEnergyLJ = calculatePotentialEnergy (particles, Time, dt, LJPotential);

    std::vector<double> totalEnergy;
    totalEnergy = calculateTotalEnergy (particles, Time, dt, LJPotential);

    return 0;
}