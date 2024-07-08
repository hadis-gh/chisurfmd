#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/Particle.h"

template<typename T>
class LennardJonesForce {
public:
    LennardJonesForce(T epsilon, T sigma) : epsilon(epsilon), sigma(sigma) {}

    T operator()(const T r) const {
        if (r == 0) return 0;
        const T sigma6 = std::pow(sigma, 6);
        const T sigma12 = sigma6 * sigma6;
        const T r6 = std::pow(r, 6);
        const T r12 = r6 * r6;
        return 48.0 * epsilon * (sigma12 / (r12 * r) - 0.5 * sigma6 / (r6 * r));
    }

private:
    T epsilon;
    T sigma;
};

template<typename T>
class LennardJonesPotential {
public:
    LennardJonesPotential(T epsilon, T sigma) : epsilon(epsilon), sigma(sigma) {}

    T operator()(const T r) const {
        if (r == 0) return 0;
        const T sigma6 = std::pow(sigma, 6);
        const T sigma12 = sigma6 * sigma6;
        const T r6 = std::pow(r, 6);
        const T r12 = r6 * r6;
        return 4.0 * epsilon * (sigma12 / r12 - 0.5 * sigma6 / r6);
    }

private:
    T epsilon;
    T sigma;
};