#pragma once

#include <iostream>
#include <cmath>
#include <tuple>
#include "Vec.h"
#include "Circle.h"

template<typename T>
T IntersectionCirclesDistance(const Circle<T> &staticCircle, const Circle<T> &movingCircle, const Vec<T> &direction)
{
    const T A = direction.abs2();
    const T B = static_cast<T>(2) * direction.dot(movingCircle.c - staticCircle.c);
    const T C = (movingCircle.c - staticCircle.c).abs2() - (staticCircle.r + movingCircle.r)*(staticCircle.r + movingCircle.r);
    
    const T discriminant = B * B - 4 * A * C;
    if (discriminant < 0) {
        return NAN;
    }

    const T t1 = (-B - std::sqrt(discriminant)) / (2 * A);
    const T t2 = (-B + std::sqrt(discriminant)) / (2 * A);

    T t;
    if (t1 >= 0 && t2 >= 0) {
        t = std::min(t1, t2);
    } else if (t1 >= 0) {
        t = t1;
    } else if (t2 >= 0) {
        t = t2;
    } else {
        return NAN;
    }

    return t;
}

template<typename T>
Vec<T> IntersectionCirclesPoint(const Circle<T> &staticCircle, const Circle<T> &movingCircle, const Vec<T> &direction){
    T distance = IntersectionCirclesDistance(staticCircle, movingCircle, direction);
    return movingCircle.c + distance *direction;
}

template<typename T>
std::tuple <Vec<T>, T> collisionProperties(const Circle<T> &staticCircle, const Circle<T> &movingCircle, const Vec<T> &direction){
    std::tuple <Vec<T>, T> result;
    return {IntersectionCirclesPoint(staticCircle, movingCircle, direction), IntersectionCirclesDistance(staticCircle, movingCircle, direction)};
}