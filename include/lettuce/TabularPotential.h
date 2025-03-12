#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/ParticleDot.h"

#include "lettuce/interpolation/interpolation.h"
#include "lettuce/interpolation/GridSelector.h"
#include "lettuce/interpolation/IrregularSelector.h"
#include "lettuce/interpolation/UnifyingYContainer.h"

template<typename T>
class TabularPotentialForce {
public:
    using value_type = T;

    TabularPotentialForce(const std::string& potentialName, const std::string& adiosInput, adios2::IO& bpIO, const std::vector<T>& params)
        : pot(potentialName, adiosInput, bpIO, params) {}

    Vec<T, 2> operator()(const T r, const T phi1, const T phi2) const {
        const auto g = interpolateGradLinear(pot.hold(), phi1, phi2, r);
        return g[2];
    }

private:
    lettuce::TabularPotential<UGrid2IrrHold<T>> pot;
};

template<typename T>
class TabularPotentialPotential {
public:
    TabularPotentialPotential(const std::string& potentialName, const std::string& adiosInput, adios2::IO& bpIO, const std::vector<T>& params)
        : pot(potentialName, adiosInput, bpIO, params) {}

    T operator()(const T r, const T phi1, const T phi2) const {
        return interpolate(InterpolateLinear(), pot.hold(), phi1, phi2, r);
    }

private:
    lettuce::TabularPotential<UGrid2IrrHold<T>> pot;
};