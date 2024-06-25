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
Circle<T> findStopPointAll (const Circle<T> &startCircle, Vec<T> &direction, const std::vector<Circle<T>> &circles){
    std::vector<std::pair<Circle<T>, T>> pairs;
    for (const auto &c: circles){
        pairs.push_back(findStopPoint(startCircle, direction, c));
    }
    std::cout << ((pairs[3]).first).c <<std::endl; //the problem is some of them are NAN and while choosing the min value they have priority to real numbers!

    auto it = std::min_element(pairs.begin(), pairs.end(),
        [](const std::pair<Circle<T>, T>& a, const std::pair<Circle<T>, T>& b) {
            if (!std::isnan(a.second) && !std::isnan(b.second)) {
                return a.second < b.second;
            }
            else if (std::isnan(a.second)) {
                return false;
            }
            else {
                return true;
            }
        });

    return it->first;
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
        movedCircle.c += direction/static_cast<T>(100.0);
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

    // Circle<float> closeCs3 = findClosestCircle3Line(startCircle.c, direction, circles);     //finding closest one with 3line method
    // Circle<float> movedCircle = moveCircle(startCircle, direction, circles, areaSize);      //moving by steps
    // Circle<float> moved2Circle = findStopPoint(startCircle, direction, closeCs3).first;     //location by solving eq

    Circle<float> finalPos = findStopPointAll (startCircle, direction, circles);            //finding closest circle and final location using solve eq for circles

    circles.push_back(startCircle);
    circles.push_back(finalPos);

    clock_t endTime = clock();
    double timeTaken = double(endTime - startTime) / CLOCKS_PER_SEC;
    std::cout << "Time taken: " << timeTaken << " seconds\n";
    std::cout << "===========================\n";
    // std::cout << "moved circle with iteration:"<< movedCircle.c <<"\n";
    // std::cout << "moved circle with solving eq:"<< moved2Circle.c <<"\n";
    std::cout << "final Position with solving eq:"<< finalPos.c <<"\n";

    writeCricles(circles.begin(), circles.begin()+circlesNumber, "circle3Lines.txt");
    writeCricles(circles.begin()+circlesNumber, circles.end(), "circle3LinesMoving.txt");

    return 0;
}