#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include "lettuce/IntersectionLineCircle.h"
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/CircleDistribution.h"

float r = 1.0;

template<typename T>
std::pair<Circle<T>, T> findClosestCircle (const Vec<T> &startPoint, const Vec<T> &direction, const std::vector<Circle<T>> &circles){
    std::vector<size_t> circlesIndex;
    std::vector<T> circlesDistance;
    T distance = 0;
    for (size_t i=0; i<circles.size() ; ++i){
        distance = calculateIntersection(startPoint, direction, circles[i]);
        if (!std::isnan(distance)){
            circlesIndex.push_back(i);
            circlesDistance.push_back(distance);
        } 
    }
    if (circlesDistance.empty()) {
        std::cout << "No intersection found!" << std::endl;
    }

    auto minDistanceIt = std::min_element(circlesDistance.begin(), circlesDistance.end());
    T minDistance = *minDistanceIt;
    size_t minIndex = std::distance(circlesDistance.begin(), minDistanceIt);

    return {circles[circlesIndex[minIndex]], minDistance};
}

template<typename T>
bool checkCollision(const Circle<T>& c1, const Circle<T>& c2) {
    return (c1.c - c2.c).abs2() <= (c1.r + c2.r)*(c1.r + c2.r);
}

template<typename T>
Circle<T> findClosestCircle3Line (const Vec<T> &startPoint, const Vec<T> &direction, const std::vector<Circle<T>> &circles){
    std::vector<std::pair<Circle<T>, T>> pairs = {findClosestCircle({{startPoint[0] - r*direction[1], startPoint[0] + r*direction[0]}}, direction, circles),
                                                  findClosestCircle({{startPoint[0] + r*direction[0], startPoint[0] + r*direction[1]}}, direction, circles),
                                                  findClosestCircle({{startPoint[0] + r*direction[1], startPoint[0] - r*direction[0]}}, direction, circles)};

    auto it = std::min_element(pairs.begin(), pairs.end(),
                              [](const std::pair<Circle<T>, T>& a,
                                 const std::pair<Circle<T>, T>& b) {
                                  return a.second < b.second;
                              });
    return it->first;
}

template<typename T>
Circle<T> moveCircle(const Circle<T> &startCircle, const Vec<T> &direction, const std::vector<Circle<T>> &circles, const T &areaSize){
    Circle<T> movedCircle = startCircle;
    bool collision = false;
    Circle<T> targetCircle = findClosestCircle3Line (movedCircle.c, direction, circles);
    while (!collision && movedCircle.c < areaSize) {
        movedCircle.c += direction/static_cast<T>(10.0);
        for (const auto& circle : circles) {
            if (checkCollision(movedCircle, targetCircle)) {
                collision = true;
                break;
            }
        }
    }
    if (!collision)
        movedCircle.c.invalidate();
    return movedCircle;
}

template<typename T>
void writeCricles(T begin, T end, const std::string& fname){
        std::ofstream output_file(fname);
        for (auto c = begin; c != end; ++c)
            output_file << c->c << ", " << c->r <<"\n";
}

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
    Circle<float> closeCs3 = findClosestCircle3Line(startCircle.c, direction, circles);
    Circle<float> movedCircle = moveCircle(startCircle, direction, circles, areaSize);

    circles.push_back(startCircle);
    circles.push_back(movedCircle);

    clock_t endTime = clock();
    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";
    std::cout << "===========================\n";
    std::cout << "Closest Circle from 3 lines: Center" << closeCs3.c << "\n";

    writeCricles(circles.begin(), circles.begin()+circlesNumber, "circle3Lines.txt");
    writeCricles(circles.begin()+circlesNumber, circles.end(), "circle3LinesMoving.txt");

    return 0;
}