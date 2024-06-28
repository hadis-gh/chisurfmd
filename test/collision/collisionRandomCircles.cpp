#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/CircleDistribution.h"
#include "lettuce/CirclesIntersectionFuncs.h"

float r = 1.0;

int main() {
    const Circle<float> startCircle = {{{0.0, 0.0}}, r};
    Vec<float, 2> direction = {{2, 3}};
    
    const float areaSize = 100;
    const int circlesNumber = 50;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<Circle<float>> circles = distCircles (circlesNumber, areaSize, r, gen);

    clock_t startTime = clock();

    direction /= static_cast<float>(sqrt(direction.abs2())); 

    Circle<float> finalPos = findStopPointAll (startCircle, direction, circles);

    circles.push_back(startCircle);
    circles.push_back(finalPos);

    clock_t endTime = clock();
    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    writeCircles(circles.begin(), circles.begin()+circlesNumber, "circle3Lines.txt");
    writeCircles(circles.begin()+circlesNumber, circles.end(), "circle3LinesMoving.txt");

    return 0;
}