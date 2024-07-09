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
    int particlesNum = 5;

    std::random_device rd;
    std::mt19937 gen(rd());

    const double sigma = 1;
    const double epsilon = 1;
    const double mass = 1;
    const double dt = 0.01;
    const double cutoff = 20;
    const double Time = 1000.;
    const double areaL = 50.0;
    const double radius = 0.005;

    InitialParticlesConfiguration  config = RANDOM_CIRCLES;
    IntegrationMethod integrationMethod = EULER; 

    LennardJonesForce<double> LJForce(epsilon, sigma, cutoff);


    clock_t startTime = clock();

    std::vector<Particle<double>> particles(particlesNum);
    particles = initialParticles (particlesNum, areaL, radius, gen, config);
    integrate(particles, dt, Time, LJForce, integrationMethod);

    clock_t endTime = clock();

    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    return 0;
}