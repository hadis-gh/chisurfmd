#include <iostream>
#include <vector>
#include <utility>
#include "lettuce/IntersectionTwoCircles.h"
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"

//Not correct (calculating distance and position of facing two circles!)

template<typename T>
void printCollisionProperties(const Circle<T>& circle, const Circle<T>& movingCircle, const Vec<T>& direction) {
    std::tuple <Vec<T>, T> result = collisionProperties(circle, movingCircle, direction);
    std::cout << std::get<0>(result) << std::endl;
    std::cout << std::get<1>(result) << std::endl;
}

template<typename T>
std::tuple<Circle<T>, Vec<T>> findClosestCircle(const std::vector<Circle<T>> &circles, const Circle<T> &movingCircle, const Vec<T> &direction){
    T minDistance = 0;
    Circle<T> closestCircle;
    Vec<T> stopPosition;

    for (Circle<T> c : circles){
        std::tuple <Vec<T>, T> collisionProperties = collisionProperties(c, movingCircle, direction);

        T distance = std::get<1>(collisionProperties);
        if (distance < minDistance) {
            minDistance = distance;
            closestCircle = c;
            stopPosition = (std::get<0>(collisionProperties)).c;
        }
    }

    return {closestCircle, stopPosition};
}


int main(){

    const Circle<float> movingCircle = {{{0.0, 0.0}}, 1.0};
    const Vec<float> direction ({1.0, 1.0});

    const Circle<float> circle1 = {{{5.0, 4.0}}, 1.0};
    const Circle<float> circle2 = {{{2.0, 1.0}}, 1.0};
    const Circle<float> circle3 = {{{20.0, 21.0}}, 1.0};

    const std::vector<Circle<float>> circles = {circle1, circle2, circle3};
/*
    printCollisionProperties(circle1, movingCircle, direction);

    std::tuple <Circle<float>, Vec<float>> result = findClosestCircle(circles, movingCircle, direction);

    std::cout << "Closest circle: " << (std::get<0>(result)).c << std::endl;
    std::cout << "Stop position: " << std::get<1>(result) << std::endl;
*/
    return 0;
}