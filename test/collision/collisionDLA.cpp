#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <random>
#include <algorithm>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/CircleDistribution.h"
#include "lettuce/CirclesIntersectionFuncs.h"


int main() {

    const float r = 1.0;
    const float width = 50.0;
    const int shootNum = 150;

    std::vector<Circle<float>> finalCircles;

    std::random_device rd;
    std::mt19937 gen(rd());

    finalCircles = distCirclesDLA (shootNum, width, r, gen);

    writeCircles(finalCircles.begin(), finalCircles.end(), "circlesDLAtest.txt");

    return 0;
}
