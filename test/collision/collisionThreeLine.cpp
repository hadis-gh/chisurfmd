#include <iostream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <stdexcept> //runtime_error
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
        throw std::runtime_error("No intersection found");
    }

    auto minDistanceIt = std::min_element(circlesDistance.begin(), circlesDistance.end());
    T minDistance = *minDistanceIt; // Dereference iterator to get value
    size_t minIndex = std::distance(circlesDistance.begin(), minDistanceIt);

    return {circles[circlesIndex[minIndex]], minDistance};
}

template<typename T>
Circle<T> findClosestCircle3Line (const Vec<T> &startPoint, const Vec<T> &direction, const std::vector<Circle<T>> &circles){
    Vec<T> startPoint2 = {{startPoint[0]+direction[0], startPoint[1]-direction[1]}};
    Vec<T> startPoint3 = {{startPoint[0]-direction[0], startPoint[1]+direction[1]}};
    std::vector<std::pair<Circle<T>, T>> pairs = {findClosestCircle (startPoint, direction, circles), 
                                                  findClosestCircle (startPoint2, direction, circles), 
                                                  findClosestCircle (startPoint3, direction, circles)};

    auto it = std::min_element(pairs.begin(), pairs.end(),
                              [](const std::pair<Circle<T>, T>& a,
                                 const std::pair<Circle<T>, T>& b) {
                                  return a.second < b.second;
                              });
    return it->first;
}


int main() {
    const Circle<float> startCircle = {{{0.0, 0.0}}, 1.0};
    const Vec<float> direction = {{1.0, 1.0}};

    const std::vector<Circle<float>> circles = {{{{5.57, 5.3}}, 1.0},
                                                {{{6.44, 5.60}}, 1.0},
                                                {{{3.34, 3.89}}, 1.0}, 
                                                {{{7.22, 7.50}}, 1.0}};
    clock_t startTime = clock();

    // Circle<float> closest = findClosestCircle(startCircle.c, direction, circles);
    Circle<float> closeCs3 = findClosestCircle3Line(startCircle.c, direction, circles);
    clock_t endTime = clock();
    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";

    std::cout << "===========================\n";
    // std::cout << "Closest Circle: " << closest.c << "\n";

    std::cout << "Closest Circle from 3 lines: Center" << closeCs3.c << "\n";
    return 0;
}
