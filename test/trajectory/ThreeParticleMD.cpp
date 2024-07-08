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

const double sigma = 1;
const double epsilon = 1;
const double mass = 1;
const double dt = 0.01;

int main(){
    Particle<double> P1{{{0.0, 0.0}}, {{0.0, 0.0}}, mass};
    Particle<double> P2{{{1.2, 0.0}}, {{0.0, 0.0}}, mass};
    Particle<double> P3{{{0.0, 1.2}}, {{0.0, 0.0}}, mass};
    Particle<double> P4{{{0.0, 5.1}}, {{0.0, 0.0}}, mass};

    std::vector<Particle<double>> particles {P1, P2, P3, P4};

    double Time = 1000.;
    int numSteps = Time/dt;

    LennardJonesForce<double> LJForce(epsilon, sigma);

    IntegrationMethod method = EULER; 

    integrate(particles, dt, Time, LJForce, method);

    return 0;
}