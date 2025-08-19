#pragma once

#include <cmath>
#include "lettuce/core/Vec.h"

template <typename T>
auto DummyForceFunc2D() {
    return [](auto&, auto&, auto& dr, T r) {
        // return Vec<T, 2>{{std::abs(dr[0]), std::abs(dr[1])}};
        // if (dr[0] < 0) dr[0] = -dr[0];
        // if (dr[1] < 0) dr[1] = -dr[1];
        return Vec<T, 2>{{dr[0], dr[1]}};
    };
}

template <typename T>
auto LennardJonesForce2D(T epsilon, T sigma) {
    return [=](auto&, auto&, const Vec<T,2>& dr, T r) {
        if (r <= std::numeric_limits<T>::epsilon()) {
            return Vec<T,2>{{0,0}}; // avoid singularity
        }

        T sr6 = std::pow(sigma / r, 6);
        T coeff = (24 * epsilon / (r * r)) * (2 * sr6 * sr6 - sr6);

        return Vec<T,2>{{coeff * dr[0], coeff * dr[1]}};
    };
}
