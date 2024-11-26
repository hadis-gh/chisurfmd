#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/ParticleDot.h"

template<typename T>
class LennardJonesForce {
public:
    using value_type = T;
    
    LennardJonesForce(T epsilon, T sigma, T cutoff) 
        : epsilon(epsilon), sigma(sigma), cutoff(cutoff), sigma6(std::pow(sigma, 6)), sigma12(sigma6 * sigma6) {}

    T operator()(const T r) const {
        if (r == 0) return 0;
        if (r > cutoff) return 0;
        
        const T r6 = std::pow(r, 6);
        const T r12 = r6 * r6;
        return 48.0 * epsilon * (sigma12 / (r12 * r) - 0.5 * sigma6 / (r6 * r));
    }

private:
    T epsilon;
    T sigma;
    T cutoff;
    T sigma6;
    T sigma12;
};

template<typename T>
class LennardJonesPotential {
public:
    LennardJonesPotential(T epsilon, T sigma, T cutoff)
        : epsilon(epsilon), sigma(sigma), cutoff(cutoff), sigma6(std::pow(sigma, 6)), sigma12(sigma6 * sigma6) {}

    T operator()(const T r) const {
        if (r == 0) return 0;
        if (r > cutoff) return 0;
        
        const T r6 = std::pow(r, 6);
        const T r12 = r6 * r6;
        return 4.0 * epsilon * (sigma12 / r12 - sigma6 / r6);
    }

private:
    T epsilon;
    T sigma;
    T cutoff;
    T sigma6;
    T sigma12;
};