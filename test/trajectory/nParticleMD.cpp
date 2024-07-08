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
// const double cutoff = 20;
double Time = 1000.;


int main(){
    int numParticle = 5;
    std::vector<Particle<double>> particles(numParticle);

    std::random_device rd;
    std::mt19937 gen(rd());

    const double areaL = 5;
    const double radius = 0.005;

    std::vector <Circle<double>> circles = distRandomCircles(numParticle, areaL, radius, gen);

    for (int i = 0; i< circles.size(); ++i){
        particles[i].r = circles[i].c; 
    }
    
    clock_t startTime = clock();

    LennardJonesForce<double> LJForce(epsilon, sigma);

    IntegrationMethod method = VELOCITY_VERLET; 

    integrate(particles, dt, Time, LJForce, method);

    clock_t endTime = clock();
    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    return 0;
}