#include <cmath>
#include <iostream>

#include "chisurfmd/core/ParticleOriented.h"
#include "chisurfmd/potential/OrientedLJ.h"


// ==============================================================================
// TEST HELPERS
// ==============================================================================

int failures = 0;

bool near(
    double actual,
    double expected,
    double tolerance = 1e-6
) {
    return std::abs(actual - expected) < tolerance;
}

void expectNear(
    double actual,
    double expected,
    const char* message,
    double tolerance = 1e-6
) {
    if (!near(actual, expected, tolerance)) {
        std::cerr
            << "FAILED: " << message
            << "\n  expected: " << expected
            << "\n  actual:   " << actual
            << '\n';

        ++failures;
    }
}


// ==============================================================================
// LJ MINIMUM
// ==============================================================================

void testLJMinimum()
{
    using Particle = ParticleOriented<double>;

    const double epsilon = 1.0;
    const double sigma = 1.0;
    const double cutoff = 3.0;

    // angularScale = 0 gives the ordinary Lennard-Jones potential.
    OrientedLJPotential<Particle> potential(
        epsilon,
        sigma,
        cutoff,
        2,
        0.0,
        0.0
    );

    OrientedLJForce<Particle> force(
        epsilon,
        sigma,
        cutoff,
        2,
        0.0,
        0.0
    );

    Particle p1;
    Particle p2;

    const double r =
        std::pow(2.0, 1.0 / 6.0) * sigma;

    Vec<double, 2> dr{{r, 0.0}};

    const double energy =
        potential(p1, p2, dr, r);

    const auto f =
        force(p1, p2, dr, r);

    expectNear(
        energy,
        -epsilon,
        "LJ energy at equilibrium distance"
    );

    expectNear(
        f[0],
        0.0,
        "LJ force at equilibrium distance"
    );
}


// ==============================================================================
// RADIAL FORCE = -dU/dr
// ==============================================================================

void testRadialForce()
{
    using Particle = ParticleOriented<double>;

    OrientedLJPotential<Particle> potential(
        1.0,    // epsilon
        1.0,    // sigma
        3.0,    // cutoff
        2,      // phi order
        0.5,    // angular scale
        0.3     // phase
    );

    OrientedLJForce<Particle> force(
        1.0,
        1.0,
        3.0,
        2,
        0.5,
        0.3
    );

    Particle p1;
    Particle p2;

    p1.phi = 0.2;
    p2.phi = 0.8;

    const double r = 1.4;
    const double h = 1e-6;

    auto energyAt = [&](double distance) {
        Vec<double, 2> dr{{distance, 0.0}};

        return potential(
            p1,
            p2,
            dr,
            distance
        );
    };

    const double dUdr =
        (
            energyAt(r + h)
            - energyAt(r - h)
        ) / (2.0 * h);

    Vec<double, 2> dr{{r, 0.0}};

    const auto f =
        force(p1, p2, dr, r);

    expectNear(
        f[0],
        -dUdr,
        "radial force agrees with -dU/dr"
    );

    expectNear(
        f[1],
        0.0,
        "radial force has no y component"
    );
}


// ==============================================================================
// TORQUE = -dU/dphi
// ==============================================================================

void testTorque()
{
    using Particle = ParticleOriented<double>;

    OrientedLJPotential<Particle> potential(
        1.0,
        1.0,
        3.0,
        2,
        0.5,
        0.3
    );

    OrientedLJForce<Particle> force(
        1.0,
        1.0,
        3.0,
        2,
        0.5,
        0.3
    );

    Particle p1;
    Particle p2;

    p1.phi = 0.2;
    p2.phi = 0.8;

    const double r = 1.4;
    const double h = 1e-6;

    Vec<double, 2> dr{{r, 0.0}};

    auto energyAtPhi2 = [&](double phi2) {
        auto shiftedParticle = p2;
        shiftedParticle.phi = phi2;

        return potential(
            p1,
            shiftedParticle,
            dr,
            r
        );
    };

    const double dUdPhi2 =
        (
            energyAtPhi2(p2.phi + h)
            - energyAtPhi2(p2.phi - h)
        ) / (2.0 * h);

    const auto f =
        force(p1, p2, dr, r);

    expectNear(
        f[2],
        -dUdPhi2,
        "torque agrees with -dU/dphi2"
    );
}


// ==============================================================================
// CUTOFF
// ==============================================================================

void testCutoff()
{
    using Particle = ParticleOriented<double>;

    OrientedLJPotential<Particle> potential(
        1.0, 1.0, 3.0, 2, 0.5, 0.0
    );

    OrientedLJForce<Particle> force(
        1.0, 1.0, 3.0, 2, 0.5, 0.0
    );

    Particle p1;
    Particle p2;

    const double r = 3.1;

    Vec<double, 2> dr{{r, 0.0}};

    expectNear(
        potential(p1, p2, dr, r),
        0.0,
        "potential is zero beyond cutoff"
    );

    const auto f =
        force(p1, p2, dr, r);

    expectNear(f[0], 0.0, "Fx is zero beyond cutoff");
    expectNear(f[1], 0.0, "Fy is zero beyond cutoff");
    expectNear(f[2], 0.0, "torque is zero beyond cutoff");
}


// ==============================================================================
// MAIN
// ==============================================================================

int main()
{
    testLJMinimum();
    testRadialForce();
    testTorque();
    testCutoff();

    if (failures > 0) {
        std::cerr
            << failures
            << " test(s) failed.\n";

        return 1;
    }

    std::cout
        << "All OrientedLJ tests passed.\n";

    return 0;
}