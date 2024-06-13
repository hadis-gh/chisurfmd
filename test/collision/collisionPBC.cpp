#include <iostream>
#include <random>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/IntersectionLineCircle.h"
#include "lettuce/CircleDistribution.h"

template <typename T>
void movePoint (Vec<float> &point, const Vec<T> &direction,const T &length){
    point += direction;

    if (point[0] < 0)       point [0] += length;
    if (point[0] >= length) point[0] -= length;
    if (point[1] < 0)       point [1] += length;
    if (point[1] >= length) point[1] -= length;
}

int main() {
    Vec<float> originPoint ({1.0, 1.0});
    Vec<float> moveDirection ({1.0, 1.0});

    float length = 50.0;
    int circleNum = 30;
    float circleR = 0.5;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector <Circle<float>> circles = distCircles(circleNum, length+ 2* circleR, circleR, gen);
    // circles.c -= 2* circles.r;

    bool foundIntersection = false;
    Vec<float> intersectionPoint;

    int maxIterations = 1000;
    int iteration = 0;

    while (!foundIntersection && iteration < maxIterations){
        for(auto c: circles){
            intersectionPoint = calculateIntersectionPoint(originPoint, moveDirection, c);

            if (!std::isnan(intersectionPoint[0])){
                foundIntersection = true;
                break;
            }
        }
        if (!foundIntersection){
            movePoint(originPoint, moveDirection, length); 
            iteration++;
        }
    }
    
    std::cout << "Intersection point: " << intersectionPoint << std::endl;
    std::cout << "Number of iteration: " << iteration << std::endl;

    return 0;
}

    // int maxIterations = 1000; // add a maximum number of iterations
    // int iteration = 0;

    // do {
    //     for (auto c : circles) {
    //         intersectionPoint = calculateIntersectionPoint(originPoint, moveDirection, c);
    //         if (!std::isnan(intersectionPoint[0]) && !std::isinf(intersectionPoint[0])) {
    //             break; // found an intersection point, exit the loop
    //         }
    //     }

    //     if (std::isnan(intersectionPoint[0]) || std::isinf(intersectionPoint[0])) {
    //         movePoint(originPoint, moveDirection, length);
    //         iteration++;
    //     }

    // } while (iteration < maxIterations && std::isnan(intersectionPoint[0]));

    // if (std::isnan(intersectionPoint[0])) {
    //     std::cerr << "No intersection point found after " << maxIterations << " iterations." << std::endl;
    // } else {
    //     std::cout << "Intersection point: " << intersectionPoint << std::endl;
    // }