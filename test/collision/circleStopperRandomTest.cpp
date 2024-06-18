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
void moveCircle(Circle<T>& movingCircle, const std::vector<Circle<T>>& circles, T angle, T stepSize, const T &areaSize) {
    T radian = angle * M_PI / 180.0;
    T dx = std::cos(radian) * stepSize;
    T dy = std::sin(radian) * stepSize;

    bool collision = false;

    while (!collision && movingCircle.c.x <areaSize && movingCircle.c.y <areaSize) {
        movingCircle.c[0] += dx;
        movingCircle.c[1] += dy;

        for (const auto& circle : circles) {
            if (checkCollision(movingCircle, circle)) {
                collision = true;
                break;
            }
        }
    }
}

int main() {
    Circle<float> movingCircle = {{{0.0, 0.0}}, 1.0};
    const float angle = 45.0;
    const float stepSize = 0.1;
    const float areaSize = 100;

    const int circlesNumber = 50;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector<Circle<float>> circles = distCircles (circlesNumber, areaSize, movingCircle.r, gen);

    clock_t startTime = clock();

    moveCircle(movingCircle, circles, angle, stepSize, areaSize);
    circles.push_back(movingCircle);

    clock_t endTime = clock();
    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    std::ofstream output_file("circleStopperRandom.txt");
    for (auto c: circles)
        output_file << c.c <<"\n";    
    return 0;
}
