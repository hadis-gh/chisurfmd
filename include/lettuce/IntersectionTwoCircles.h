#pragma once

#include <iostream>
#include <cmath>
#include "Vec.h"
#include "Circle.h"

template<typename T>
Vec<T> findIntersectionCircles(const Circle<T> &staticCircle, const Circle<T> &movingCircle, const Vec<T> &direction)
{
    Vec<T> intersection;
    using a = staticCircle.c.x;
    using b = staticCircle.c.y;
    using r = staticCircle.r;

    using x0 = movingCircle.c.x;
    using y0 = movingCircle.c.y;
    using R = movingCircle.r;
    using dx = direction.x;
    using dy = direction.y;

    const T A = direction.abs2();
    const T B = 2 * (dx * (x0 - a) + dy * (y0 - b));
    // const T B = 2 * direction * (movingCircle.c - staticCircle.c);

    const T C = (x0 - a) * (x0 - a) + (y0 - b) * (y0 - b) - (R + r) * (R + r);
    //const T C = movingCircle.c * staticCircle.c - (R + r) * (R + r);
    

    const T discriminant = B * B - 4 * A * C;
    if (discriminant < 0) {
        return NAN;
    }

    const T t1 = (-B - std::sqrt(discriminant)) / (2 * A);
    const T t2 = (-B + std::sqrt(discriminant)) / (2 * A);

    const T t;
    if (t1 >= 0 && t2 >= 0) {
        t = std::min(t1, t2);
    } else if (t1 >= 0) {
        t = t1;
    } else if (t2 >= 0) {
        t = t2;
    } else {
        return NAN;
    }

    intersection.x = movingCircle.c + t *direction;

    return intersection;
}