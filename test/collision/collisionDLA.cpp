#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <random>
#include <algorithm>
#include "lettuce/ParticleDot.h"
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/CircleDistribution.h"
#include "lettuce/CirclesIntersectionFuncs.h"


int main() {

    const float r = 1.0;
    const float width = 50.0;
    const int shootNum = 150;

    std::vector<Circle<float>> finalCircles;
    std::vector<ParticleDot<float>> finalParticles;

    std::random_device rd;
    std::mt19937 gen(rd());

    finalCircles = distCirclesDLA (shootNum, width, r, gen);
    // finalParticles = distParticleDLA (shootNum, width, r, gen);
    finalParticles = distParticleDLA<ParticleDot<float>>(shootNum, width, r, gen);
    writeCircles(finalCircles.begin(), finalCircles.end(), "circlesDLAtest.txt");
    writeParticle(finalParticles, r, "particlesDLAtest.txt");

    return 0;
}
