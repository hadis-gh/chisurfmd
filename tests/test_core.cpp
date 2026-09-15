#include <cmath>
#include <iostream>

#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/ParticleDot.h"
#include "chisurfmd/core/ParticleOriented.h"


// ==============================================================================
// TEST HELPERS
// ==============================================================================

bool near(double a, double b, double tolerance = 1e-12)
{
    return std::abs(a - b) < tolerance;
}

int failures = 0;

void expect(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << "FAILED: " << message << '\n';
        ++failures;
    }
}


// ==============================================================================
// VEC TESTS
// ==============================================================================

void testVec()
{
    Vec<double, 2> a{{1.0, 2.0}};
    Vec<double, 2> b{{3.0, 4.0}};

    const auto sum = a + b;

    expect(
        near(sum[0], 4.0) &&
        near(sum[1], 6.0),
        "Vec addition"
    );

    const auto difference = b - a;

    expect(
        near(difference[0], 2.0) &&
        near(difference[1], 2.0),
        "Vec subtraction"
    );

    expect(
        near(a.abs2(), 5.0),
        "Vec squared magnitude"
    );
}


// ==============================================================================
// PARTICLE TESTS
// ==============================================================================

void testDegreesOfFreedom()
{
    expect(
        degreesOfFreedom<ParticleDot<double>>() == 2,
        "ParticleDot has 2 degrees of freedom"
    );

    expect(
        degreesOfFreedom<ParticleOriented<double>>() == 3,
        "ParticleOriented has 3 degrees of freedom"
    );
}


// ==============================================================================
// PERIODIC BOUNDARY TEST
// ==============================================================================

void testPeriodicBoundary()
{
    ParticleDot<double> particle;

    particle.position[0] = 10.2;
    particle.position[1] = -0.3;

    implementPBC(particle, 10.0);

    expect(
        near(particle.position[0], 0.2),
        "PBC wraps x coordinate"
    );

    expect(
        near(particle.position[1], 9.7),
        "PBC wraps y coordinate"
    );
}


// ==============================================================================
// MAIN
// ==============================================================================

int main()
{
    testVec();
    testDegreesOfFreedom();
    testPeriodicBoundary();

    if (failures > 0) {
        std::cerr
            << failures
            << " test(s) failed.\n";

        return 1;
    }

    std::cout << "All core tests passed.\n";

    return 0;
}