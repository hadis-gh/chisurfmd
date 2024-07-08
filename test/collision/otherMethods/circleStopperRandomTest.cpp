#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <ctime>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/CircleDistribution.h"

//with random circles on the surface - but result needs work/ circle distribution does not seem correct

template<typename T>
bool checkCollision(const Circle<T>& c1, const Circle<T>& c2) {
    return (c1.c - c2.c).abs2() <= (c1.r + c2.r)*(c1.r + c2.r);
}

template<typename T>
auto moveCircle(Circle<T> movingCircle, const std::vector<Circle<T>>& circles, T angle, T stepSize, const T &areaSize) {
    T radian = angle * M_PI / 180.0;
    T dx = std::cos(radian) * stepSize;
    T dy = std::sin(radian) * stepSize;

    bool collision = false;

    while (!collision && movingCircle.c < areaSize) {
        movingCircle.c[0] += dx;
        movingCircle.c[1] += dy;

        for (const auto& circle : circles) {
            if (checkCollision(movingCircle, circle)) {
                collision = true;
                break;
            }
        }
    }
    if (!collision)
        movingCircle.c.invalidate();
    return movingCircle;
}

int main() {
    const float angle = 45.0;
    const float stepSize = 0.1;
    const float radius = 1.;
    const float areaSize = 100;

    const int circlesNumber = 50;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<Circle<float>> circles = distRandomCircles (circlesNumber, areaSize, radius, gen);

    clock_t startTime = clock();

    std::uniform_real_distribution<> dis(0, areaSize);

    for (int i = 0; i < 10; ++i)
    {
        auto movingCircle = placeRandomCircle(circles, dis, radius, gen);
        if(std::isnan(movingCircle.c[0]))
            continue;
        movingCircle = moveCircle(movingCircle, circles, 45.f, stepSize, areaSize);
        if(std::isnan(movingCircle.c[0]))
            continue;
        circles.push_back(movingCircle);
    }

    auto movingCircle = circles [0];
    movingCircle.c[0] += movingCircle.r / 2;
    moveCircle(movingCircle, circles, angle, stepSize, areaSize);
    circles.push_back(movingCircle);

    clock_t endTime = clock();
    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    writeCircles(circles.begin(), circles.begin()+circlesNumber, "circleStopperRandom.txt");

    writeCircles(circles.begin()+circlesNumber, circles.end(), "circleStopperRandomMoving.txt");
    
    return 0;
}
