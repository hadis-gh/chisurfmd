#include <iostream>
#include <fstream>
#include <random>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/IntersectionLineCircle.h"
#include "lettuce/CircleDistribution.h"

template <typename T>
void movePoint (Vec<float> &point, const Vec<T> &direction,const T &length){
    point += direction;

    if (point[0] < 0)       point[0] += length;
    if (point[0] >= length) point[0] -= length;
    if (point[1] < 0)       point[1] += length;
    if (point[1] >= length) point[1] -= length;
}


int main() {
    Vec<float> originPoint ({20.0, 20.0});
    Vec<float> startPoint = originPoint;
    Vec<float> moveDirection ({1.0, 3.0});

    float length = 50.0;
    int circleNum = 25;
    float circleR = 1;

    std::random_device rd;
    std::mt19937 gen(rd());

    std::vector <Circle<float>> circles = distCirclesPBC(circleNum, length, circleR, gen);

    bool foundIntersection = false;
    Vec<float> intersectionPoint;

    int maxIterations = 1000;
    int iteration = 0;

    while (!foundIntersection && iteration < maxIterations){
        for(auto c: circles){
            intersectionPoint = calculateIntersectionPoint(startPoint, moveDirection, c);

            if (!std::isnan(intersectionPoint[0])){
                foundIntersection = true;
                break;
            }
        }
        if (!foundIntersection){
            movePoint(startPoint, moveDirection, length); 
            iteration++;
        }
    }
    
    std::cout << "Intersection point: " << intersectionPoint << std::endl;
    std::cout << "Number of iteration: " << iteration << std::endl;

    std::ofstream outputFile("collisionPBC_info.txt");
    outputFile << "originPoint= "<< originPoint
               << "\nstartPoint= "<< startPoint
               << "\nmovementDirection= "<< moveDirection
               << "\nlength= "<< length
               << "\ncirclesRadius= "<< circleR

               << "\n\nintersection= "<< intersectionPoint
               << "\nNumberOfIteration= " << iteration;

    outputFile << "\n\ncirclesCenters= [";
    for (auto c: circles){
        outputFile << c.c <<", ";
    }
    outputFile << "]"<< std::endl;

    return 0;
}