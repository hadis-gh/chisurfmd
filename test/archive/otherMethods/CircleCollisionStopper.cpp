#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <ctime>
#include "lettuce/core/Vec.h"
#include "lettuce/core/Circle.h"

//correct one for test the program

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

    while (!collision && movingCircle.c[0] <areaSize && movingCircle.c[1]<areaSize) {
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
    const float areaSize = 20;

    std::vector<Circle<float>> circles = {{{{5.57, 11.3}}, 1.0},
                                          {{{6.44, 5.60}}, 1.0},
                                          {{{3.34, 4.89}}, 1.0}, 
                                          {{{7.22, 7.50}}, 1.0}};
    clock_t startTime = clock();

    moveCircle(movingCircle, circles, angle, stepSize, areaSize);
    circles.push_back(movingCircle);

    clock_t endTime = clock();
    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    std::cout << "===========================" << "\n";
    std::cout << "Circles: " << "\n";
    for (const auto& c: circles)
        std::cout << c.c <<"\n";
    return 0;
}
