#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include "lettuce/IntersectionLineCircle.h"
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"

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

// template<typename T>
// Circle<T> findClosestCircle3Line (const Vec<T> &startPoint, const Vec<T> &direction, const std::vector<Circle<T>> &circles){
//     std::vector<std::pair<Circle<T>, T>> pairs = {findClosestCircle({{startCircle.c[0] + r*dx/h, startCircle.c[0] + r*dy/h}}, direction, circles),
//                                                   findClosestCircle({{startCircle.c[0] + r*dy/h, startCircle.c[0] - r*dx/h}}, direction, circles),
//                                                   findClosestCircle({{startCircle.c[0] - r*dy/h, startCircle.c[0] + r*dx/h}}, direction, circles)};

//     auto it = std::min_element(pairs.begin(), pairs.end(),
//                               [](const std::pair<Circle<T>, T>& a,
//                                  const std::pair<Circle<T>, T>& b) {
//                                   return a.second < b.second;
//                               });
//     return it->first;
// }

float r = 1.0;

int main() {
    const Circle<float> startCircle = {{{-1.0, 0.0}}, r};
    const Vec<float> direction = {{1.0, 1.0}};

    const auto h = static_cast<float>(sqrt(direction.abs2()));    
    const auto dx = direction[0];
    const auto dy = direction[1];

    const std::vector<Circle<float>> circles = {{{{6.52, 7.26}}, r},
                                                {{{2.00, 6.00}}, r},
                                                {{{1.63, 3.61}}, r}, 
                                                {{{4.59, 2.26}}, r}};
    clock_t startTime = clock();

    // Circle<float> closeCs3 = findClosestCircle3Line(startCircle.c, direction, circles);
    
    std::pair<Circle<float>, float> closeTestLine1 = findClosestCircle({{startCircle.c[0] + r*dx/h, startCircle.c[0] + r*dy/h}}, direction, circles);
    std::pair<Circle<float>, float> closeTestLine2 = findClosestCircle({{startCircle.c[0] + r*dy/h, startCircle.c[0] - r*dx/h}}, direction, circles);
    std::pair<Circle<float>, float> closeTestLine3 = findClosestCircle({{startCircle.c[0] - r*dy/h, startCircle.c[0] + r*dx/h}}, direction, circles);

    clock_t endTime = clock();
    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";
    std::cout << "===========================\n";
    // std::cout << "Closest Circle from 3 lines: Center" << closeCs3.c << "\n";



    std::cout << "Closest Circle from 1 lines: Center" << closeTestLine1.first.c << "\n";
    std::cout << "Closest Circle from 1 lines: " << closeTestLine1.second << "\n";

    std::cout << "Closest Circle from 1 lines: Center" << closeTestLine2.first.c << "\n";
    std::cout << "Closest Circle from 1 lines: " << closeTestLine2.second << "\n";

    std::cout << "Closest Circle from 1 lines: Center" << closeTestLine3.first.c << "\n";
    std::cout << "Closest Circle from 1 lines: " << closeTestLine3.second << "\n";
    return 0;
}