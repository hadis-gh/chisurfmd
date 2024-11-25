#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/ParticleDot.h"

template<typename T>
class LennardJonesOrientedForce {
public:
    using value_type = T;
    
    LennardJonesOrientedForce(T epsilon, T sigma, T cutoff, T phiConst) 
        : epsilon(epsilon), sigma(sigma), cutoff(cutoff), sigma6(std::pow(sigma, 6)), sigma12(sigma6 * sigma6), phiConst(phiConst) {}

    Vec<T, 2> operator()(const T r, const T deltaPhi) const {
        if (r == 0 || r > cutoff) return {{0, 0}};
        T phi_c = 2.0;

        const T radialForce = 4.0 * epsilon * phiConst * cos(phi_c * deltaPhi) * (
            -12.0 * sigma12 / std::pow(r, 13) + 3.0 * sigma6 / std::pow(r, 7)
        );

        const T angularForce = 4.0 * epsilon * phiConst * (sigma12 / std::pow(r, 12) - sigma6 / (2.0 * std::pow(r, 6))) * phi_c * sin(phi_c * deltaPhi);

        return {{radialForce, angularForce}};
    }
    
private:
    T epsilon;
    T sigma;
    T cutoff;
    T sigma6;
    T sigma12;
    T phiConst;
};

template<typename T>
class LennardJonesOrientedPotential {
public:
    LennardJonesOrientedPotential(T epsilon, T sigma, T cutoff, T phiConst)
        : epsilon(epsilon), sigma(sigma), cutoff(cutoff), sigma6(std::pow(sigma, 6)), sigma12(sigma6 * sigma6), phiConst(phiConst) {}

    T operator()(const T r, const T deltaPhi) const {
        if (r == 0) return 0;
        if (r > cutoff) return 0;
        
        const T r6 = std::pow(r, 6);
        const T r12 = r6 * r6;
        return 4.0 * epsilon * (sigma12 / r12 - 0.5 * sigma6 / r6) * phiConst * cos(4 * deltaPhi);
    }

private:
    T epsilon;
    T sigma;
    T cutoff;
    T sigma6;
    T sigma12;
    T phiConst;
};