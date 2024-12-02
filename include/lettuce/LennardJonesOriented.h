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
        : epsilon(epsilon), sigma(sigma), cutoff(cutoff), sigma6(std::pow(sigma, 6))
            , sigma12(sigma6 * sigma6), phiConst(phiConst), rotationalOrder(rotationalOrder) {}

    Vec<T, 2> operator()(const T r, const T deltaPhi) const {
        if (r == 0 || r > cutoff) return {{0, 0}};
        // T phi_c = 2.0;

        // const T radialForce = 4.0 * epsilon * phiConst * cos(phi_c * deltaPhi) * (
        //     -12.0 * sigma12 / std::pow(r, 13) + 3.0 * sigma6 / std::pow(r, 7)
        // );

        // const T angularForce = 4.0 * epsilon * phiConst * (sigma12 / std::pow(r, 12) - sigma6 / (2.0 * std::pow(r, 6))) * phi_c * sin(phi_c * deltaPhi);

        // return {{radialForce, angularForce}};

        const T rotationalOrder = 2;
        const T A = 5 * std::pow(r, -6);
        const T dA_dr = - A * 6 / r;

        const T radialForce = -4.0 * epsilon * (-12.0 * std::pow(sigma / r, 12) / r
            + 6.0 * std::pow(sigma / r, 6) / r) + dA_dr * std::cos(rotationalOrder * deltaPhi);

        const T angularForce = rotationalOrder * A * std::sin(rotationalOrder * deltaPhi);

        return {{radialForce, angularForce}};
    }
    
private:
    T epsilon;
    T sigma;
    T cutoff;
    T sigma6;
    T sigma12;
    T phiConst;
    int rotationalOrder; /// set it in main
};

//***************************/ new section that should be discussed/***************************

template<typename T>
class AdditivePotentialForce {
public:
    AdditivePotentialForce(T epsilon, T sigma, T cutoff, T n) 
        : epsilon(epsilon), sigma(sigma), cutoff(cutoff), n(n) {}

    Vec<T, 2> operator()(const T r, const T deltaPhi) const {
        if (r == 0 || r > cutoff) return {{0, 0}};

        const T dA_dr = 0; // Ar

        // LennardJonesForce<T> LjForce(epsilon, sigma, cutoff);
        // const T radialForce = LjForce(r) + dA_dr * std::cos(n * deltaPhi);

        const T radialForce = -4.0 * epsilon * (-12.0 * std::pow(sigma / r, 12) / r + 6.0 * std::pow(sigma / r, 6) / r) - 30 * std::pow(r, -7) * std::cos(n * deltaPhi);

        const T angularForce = n * 5 * std::pow(r, -6) * std::sin(n * deltaPhi);

        return {{radialForce, angularForce}};
    }

private:
    T epsilon;
    T sigma;
    T cutoff;
    T n;
};

template<typename T>
class MultiplicativePotentialForce {
public:
    MultiplicativePotentialForce(T epsilon, T sigma, T cutoff, T n) 
        : epsilon(epsilon), sigma(sigma), cutoff(cutoff), n(n) {}

    Vec<T, 2> operator()(const T r, const T deltaPhi) const {
        if (r == 0 || r > cutoff) return {{0, 0}};

        // LennardJonesPotential<T> LjPotential(epsilon, sigma, cutoff);
        // const T radialPotential = LjPotential(r);

        const T radialPotential = 4.0 * epsilon * (std::pow(sigma / r, 12) - std::pow(sigma / r, 6));

        const T dRadialPotential_dr = -4.0 * epsilon * (-12.0 * std::pow(sigma / r, 12) / r + 6.0 * std::pow(sigma / r, 6) / r);
        const T radialForce = dRadialPotential_dr * std::cos(n * deltaPhi);

        const T angularForce = n * radialPotential * std::sin(n * deltaPhi);

        return {{radialForce, angularForce}};
    }

private:
    T epsilon;
    T sigma;
    T cutoff;
    T n;       // Angular periodicity
};

//******************************************************

template<typename T>
class LennardJonesOrientedPotential {
public:
    LennardJonesOrientedPotential(T epsilon, T sigma, T cutoff, T phiConst)
        : epsilon(epsilon), sigma(sigma), cutoff(cutoff), sigma6(std::pow(sigma, 6))
        , sigma12(sigma6 * sigma6), phiConst(phiConst), rotationalOrder(rotationalOrder) {}

    T operator()(const T r, const T deltaPhi) const {
        if (r == 0) return 0;
        if (r > cutoff) return 0;
        
        const T r6 = std::pow(r, 6);
        const T r12 = r6 * r6;
        return 4.0 * epsilon * (sigma12 / r12 - sigma6 / r6) * phiConst * cos(4 * deltaPhi);

        // LennardJonesPotential<T> LjPotential(epsilon, sigma, cutoff);
        // return LjPotential(r) * cos(4 * deltaPhi);
    }

private:
    T epsilon; //Member variables should start with m_ like m_epsilon
    T sigma;
    T cutoff;
    T sigma6;
    T sigma12;
    T phiConst;
    int rotationalOrder;
};