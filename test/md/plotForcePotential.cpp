#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include "lettuce/core/Vec.h"
#include "lettuce/core/Circle.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/Utilities.h"
#include "lettuce/md/Integration.h"
#include "lettuce/potential/IsotropicLJ.h"
#include "lettuce/core/CircleDistribution.h"
#include "lettuce/core/CirclesIntersectionFuncs.h"


const double sigma = 1;
const double epsilon = 1;
const double mass = 1;
const double cutoff = 20;
const double dt = 0.01;

// template<typename T>
// void writePlotData(const T &startR, const T &endR, const int &numSpace){
//     std::vector<double> distances = linspace(startR, endR, numSpace);
//     std::vector<double> ljForceValue(numSpace);
//     std::vector<double> ljPotentialValue(numSpace);

//     LennardJonesForce<ParticleDot<double>> LJForce(epsilon, sigma, cutoff);
//     LennardJonesPotential<ParticleDot<double>> LJPotential(epsilon, sigma, cutoff);

//     for (int i = 0; i < distances.size(); ++i){
//         ljForceValue[i] = LJForce(distances[i], 0);
//         ljPotentialValue[i] = LJPotential(distances[i], 0);
//     }

//     std::ofstream ljForceFile ("plotLJforce.txt");
//     for (int i= 0; i < distances.size(); ++i){
//         ljForceFile << distances[i] << " " << ljForceValue[i] << std::endl;
//     }

//     std::ofstream ljPotentialFile ("plotLJpotential.txt");
//     for (int i= 0; i < distances.size(); ++i){
//         ljPotentialFile << distances[i] << " " << ljPotentialValue[i] << std::endl;
//     }
// }


int main(){

    const double startR = 0.9 * sigma;
    const double endR = 3.0 * sigma;
    const int numSpace = 400;

    // writePlotData(startR, endR, numSpace);

    return 0;
}
