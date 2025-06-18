#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include "Vec.h"
#include "Circle.h"

template<typename T>
std::pair<Circle<T>, T> findStopPoint (const Circle<T> &startCircle, Vec<T> &direction, const Circle<T> &closestCircle){

    const T a = direction.abs2();
    const T b = static_cast<T>(2.0) * direction*(startCircle.c - closestCircle.c);
    const T c = (startCircle.c - closestCircle.c).abs2() - 4 * startCircle.r * startCircle.r;

    const T discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {{NAN, NAN}, NAN};
    }

    const T t1 = (-b + std::sqrt(discriminant)) / (2 * a);
    const T t2 = (-b - std::sqrt(discriminant)) / (2 * a);

    T t;
    if (t1 >= 0 && t2 >= 0) {
        t = std::min(t1, t2);
    } else if (t1 >= 0) {
        t = t1;
    } else if (t2 >= 0) {
        t = t2;
    } else {
        return {{NAN, NAN}, NAN};
    }
    return {{startCircle.c + t * direction}, t};
}

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
// bool checkCollision(const Circle<T>& c1, const Circle<T>& c2) {
//     return (c1.c - c2.c).abs2() <= (c1.r + c2.r)*(c1.r + c2.r);
// }

// template<typename T>
// Circle<T> findClosestCircle3Line (const Circle<T> &startPoint, const Vec<T> &direction, const std::vector<Circle<T>> &circles){
//     std::vector<std::pair<Circle<T>, T>> pairs = {findClosestCircle({{startPoint.c[0] - startPoint.r*direction[1], startPoint.c[0] + startPoint.r*direction[0]}}, direction, circles),
//                                                   findClosestCircle({{startPoint.c[0] + startPoint.r*direction[0], startPoint.c[0] + startPoint.r*direction[1]}}, direction, circles),
//                                                   findClosestCircle({{startPoint.c[0] + startPoint.r*direction[1], startPoint.c[0] - startPoint.r*direction[0]}}, direction, circles)};

//     auto it = std::min_element(pairs.begin(), pairs.end(),
//                               [](const std::pair<Circle<T>, T>& a,
//                                  const std::pair<Circle<T>, T>& b) {
//                                   return a.second < b.second;
//                               });
//     return it->first;
// }

// template<typename T>
// Circle<T> moveCircle(const Circle<T> &startCircle, const Vec<T> &direction, const std::vector<Circle<T>> &circles, const T &areaSize){
//     Circle<T> movedCircle = startCircle;
//     bool collision = false;
//     Circle<T> targetCircle = findClosestCircle3Line (movedCircle.c, direction, circles);
//     while (!collision && movedCircle.c < areaSize) {
//         movedCircle.c += direction/static_cast<T>(100.0);
//         for (const auto& circle : circles) {
//             if (checkCollision(movedCircle, targetCircle)) {
//                 collision = true;
//                 break;
//             }
//         }
//     }
//     if (!collision)
//         movedCircle.c.invalidate();
//     return movedCircle;
// }