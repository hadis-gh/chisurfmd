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

template<typename T>
void writePlotData(const T &startR, const T &endR, const int &numSpace){
    std::vector<double> distances = linspace(startR, endR, numSpace);
    std::vector<double> ljForceValue(numSpace);
    std::vector<double> ljPotentialValue(numSpace);

    LennardJonesForce<double> LJForce(epsilon, sigma);
    LennardJonesPotential<double> LJPotential(epsilon, sigma);

    for (int i = 0; i < distances.size(); ++i){
        ljForceValue[i] = LJForce(distances[i]);
        ljPotentialValue[i] = LJPotential(distances[i]);
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

    const double startR = 0.9 * sigma;
    const double endR = 3.0 * sigma;
    const int numSpace = 400;

    writePlotData(startR, endR, numSpace);

    return 0;
}
