#pragma once

#include <iostream>
#include <cmath>
#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/Circle.h"

//old functions - probably not useful anymore

template <typename T>
T calculateIntersection(const Vec<T, 2> &point,const Vec<T> &direction,const Circle<T> &circle) 
{   
    const T a = direction.dot(direction);
    const auto dist = point - circle.c;
    const T b = 2 * direction.dot(dist);
    // T b = 2 * (dx * (point[0] - circle.x) + dy * (py - circle.y));

    // const T c = (px - circle.x) * (px - circle.x) + (py - circle.y) * (py - circle.y) - circle.r * circle.r;
    const T c = dist.dot(dist) - circle.r * circle.r;
    const T discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return NAN;
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
        return NAN;
    }

    return t;
}

template <typename T>
Vec<T, 2> calculateIntersectionPoint(const Vec<T, 2> &point,const Vec<T> &direction,const Circle<T> &circle) 
{
    const auto t = calculateIntersection(point, direction, circle);
    return point + t*direction;
}

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