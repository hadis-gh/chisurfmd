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

const float r = 1.0;
const float width = 50.0;

int main() {

    Circle<float> target = {{{width / 2, width / 2}}, r};
    std::vector<Circle<float>> finalCircles;
    finalCircles.push_back(target);

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int j = 0; j < 50; ++j) {
        Circle<float> newCircle = startCircleRandom(r, width, gen);
        Vec<float> direction = shootToCenter(newCircle, width);

        Circle<float> endPoint = findStopPointAll(newCircle, direction, finalCircles);
        if (!std::isnan(endPoint.c[0])){
            finalCircles.push_back(endPoint);
        }
    }
    writeCircles(finalCircles.begin(), finalCircles.end(), "circlesDLAtest.txt");

    return 0;
}
