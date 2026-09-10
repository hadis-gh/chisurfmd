#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "chisurfmd/core/ParticleOriented.h"
#include "chisurfmd/core/Species.h"
#include "chisurfmd/md/Integration.h"
#include "chisurfmd/potential/OrientedLJ.h"


// ==============================================================================
// TEST HELPERS
// ==============================================================================

int failures = 0;

void expectLess(
    double actual,
    double limit,
    const char* message
) {
    if (actual >= limit) {
        std::cerr
            << "FAILED: " << message
            << "\n  value: " << actual
            << "\n  limit: " << limit
            << '\n';

        ++failures;
    }
}


// ==============================================================================
// ENERGY
// ==============================================================================

template<typename Particle, typename Potential>
double pairPotentialEnergy(
    const Particle& p1,
    const Particle& p2,
    double boxSize,
    const Potential& potential
) {
    Vec<double, 2> dr =
        p2.position - p1.position;

    // Minimum-image convention.
    for (int i = 0; i < 2; ++i) {
        if (dr[i] > boxSize / 2.0)
            dr[i] -= boxSize;

        else if (dr[i] < -boxSize / 2.0)
            dr[i] += boxSize;
    }

    const double r = dr.abs();

    return potential(
        p1,
        p2,
        dr,
        r
    );
}


template<typename Particle>
double kineticEnergy(
    const std::vector<Particle>& particles,
    const std::vector<Species<double>>& species
) {
    double energy = 0.0;

    for (const auto& p : particles) {

        const auto velocity =
            getGeneralizedVelocities(p);

        const double mass =
            species[p.species].mass;

        const double momentI =
            species[p.species].momentOfInertia;

        energy +=
            0.5 * mass *
            (
                velocity[0] * velocity[0]
                + velocity[1] * velocity[1]
            );

        if (velocity.size() > 2) {
            energy +=
                0.5 * momentI *
                velocity[2] * velocity[2];
        }
    }

    return energy;
}


// ==============================================================================
// MOMENTUM
// ==============================================================================

template<typename Particle>
Vec<double, 2> totalMomentum(
    const std::vector<Particle>& particles,
    const std::vector<Species<double>>& species
) {
    Vec<double, 2> momentum{{0.0, 0.0}};

    for (const auto& p : particles) {

        const double mass =
            species[p.species].mass;

        momentum += p.velocity * mass;
    }

    return momentum;
}


// ==============================================================================
// VELOCITY VERLET: NVE TEST
// ==============================================================================

void testVelocityVerletNVE()
{
    using Particle =
        ParticleOriented<double>;

    // --------------------------------------------------------------------------
    // System
    // --------------------------------------------------------------------------

    const double boxSize = 20.0;
    const double dt = 1e-3;
    const int steps = 5000;

    std::vector<Species<double>> species{
        Species<double>{
            1.0,   // mass
            1.0,   // moment of inertia
            0.5    // radius
        }
    };


    // --------------------------------------------------------------------------
    // Potential
    // --------------------------------------------------------------------------

    // angularScale = 0 reduces OrientedLJ to ordinary Lennard-Jones.
    OrientedLJPotential<Particle> potential(
        1.0,    // epsilon
        1.0,    // sigma
        3.0,    // cutoff
        2,      // phi order
        0.0,    // angular scale
        0.0     // phase
    );

    OrientedLJForce<Particle> force(
        1.0,
        1.0,
        3.0,
        2,
        0.0,
        0.0
    );


    // --------------------------------------------------------------------------
    // Initial configuration
    // --------------------------------------------------------------------------

    Particle p1;
    Particle p2;

    p1.position[0] = 9.35;
    p1.position[1] = 10.0;

    p2.position[0] = 10.65;
    p2.position[1] = 10.0;

    p1.velocity[0] = 0.0;
    p1.velocity[1] = 0.0;

    p2.velocity[0] = 0.0;
    p2.velocity[1] = 0.0;

    p1.phi = 0.0;
    p2.phi = 0.0;

    p1.omega = 0.0;
    p2.omega = 0.0;

    std::vector<Particle> particles{
        p1,
        p2
    };


    // --------------------------------------------------------------------------
    // Initial energy and momentum
    // --------------------------------------------------------------------------

    const double initialEnergy =
        kineticEnergy(particles, species)
        + pairPotentialEnergy(
            particles[0],
            particles[1],
            boxSize,
            potential
        );

    const auto initialMomentum =
        totalMomentum(particles, species);

    double maxEnergyDrift = 0.0;


    // --------------------------------------------------------------------------
    // Integration
    // --------------------------------------------------------------------------

    for (int step = 0; step < steps; ++step) {

        VelocityVerletStep(
            particles,
            species,
            dt,
            boxSize,
            force
        );

        const double energy =
            kineticEnergy(particles, species)
            + pairPotentialEnergy(
                particles[0],
                particles[1],
                boxSize,
                potential
            );

        maxEnergyDrift =
            std::max(
                maxEnergyDrift,
                std::abs(energy - initialEnergy)
            );
    }


    // --------------------------------------------------------------------------
    // Validation
    // --------------------------------------------------------------------------

    const double relativeEnergyDrift =
        maxEnergyDrift /
        std::abs(initialEnergy);

    expectLess(
        relativeEnergyDrift,
        1e-4,
        "Velocity Verlet conserves total energy"
    );

    const auto finalMomentum =
        totalMomentum(particles, species);

    expectLess(
        std::abs(
            finalMomentum[0]
            - initialMomentum[0]
        ),
        1e-10,
        "x momentum is conserved"
    );

    expectLess(
        std::abs(
            finalMomentum[1]
            - initialMomentum[1]
        ),
        1e-10,
        "y momentum is conserved"
    );

    std::cout
        << "Maximum relative energy drift: "
        << relativeEnergyDrift
        << '\n';
}


// ==============================================================================
// MAIN
// ==============================================================================

int main()
{
    testVelocityVerletNVE();

    if (failures > 0) {
        std::cerr
            << failures
            << " test(s) failed.\n";

        return 1;
    }

    std::cout
        << "All integration tests passed.\n";

    return 0;
}