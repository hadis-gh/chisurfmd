#include <iostream>
#include <fstream>
#include "lettuce/CircleDistribution.h"
#include "lettuce/ParticleDot.h"
#include "lettuce/ParticleOriented.h"

int main() {
    int particlesNum = 49;
    double areaL = 20.0;
    double radius = .8;
    double mass = 1;
    Species<double> species1 {mass, 0, radius};
    Species<double> species2 {2.0f * mass, 0, 0.5f * radius};
    std::vector<Species<double>> allSpecies {species1, species2};
    int speciesInd = 0;

    std::random_device rd;
    std::mt19937 gen(rd());

    auto particles1 = initialParticles<ParticleDot<double>>(particlesNum, allSpecies, speciesInd, areaL, gen, "RANDOM", "RRUU");
    auto particles2 = initialParticles<ParticleDot<double>>(particlesNum, allSpecies, speciesInd, areaL, gen, "RANDOM2", "RRUU");
    auto particles3 = initialParticles<ParticleDot<double>>(particlesNum, allSpecies, speciesInd, areaL, gen, "DLA", "RRUU");

    std::ofstream outputFile1("particleGeneralizedTest_RANDOM1.dat");
    std::ofstream outputFile2("particleGeneralizedTest_RANDOM2.dat");
    std::ofstream outputFile3("particleGeneralizedTest_DLA.dat");

    for (int i = 0; i < particles1.size(); ++i){
        outputFile1 << particles1[i].r[0] << " " << particles1[i].r[1] <<"\n";
        outputFile2 << particles2[i].r[0] << " " << particles2[i].r[1] <<"\n";
        outputFile3 << particles3[i].r[0] << " " << particles3[i].r[1] <<"\n";
    }

    return 0;
}