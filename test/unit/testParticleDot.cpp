#define CATCH_CONFIG_MAIN
#include <random>
#include "./include/catch_amalgamated.hpp"
#include "lettuce/core/Vec.h"
#include "lettuce/core/Species.h"
#include "lettuce/core/ParticleDot.h"

using Catch::Approx;
using Real = double;
using ParticleT = ParticleDot<Real>;
const Real areaL = 10;
const Real mass = 1;
const Real momentI = 1;
const Real radius = 0.5;
std::vector<Species<Real>> species{{mass, momentI, radius}}; // mass=1

// dr as force
auto dummyForce = [](auto&, auto&, auto& dr, double r) {
    return Vec<double,2>{{std::abs(dr[0]), std::abs(dr[1])}};
};

TEST_CASE("Default constructor") {
    ParticleT p;
    REQUIRE(p.species == 0u);
    REQUIRE(p.r[0] == 0.0);
    REQUIRE(p.r[1] == 0.0);
    REQUIRE(p.v[0] == 0.0);
    REQUIRE(p.h == 1);
    REQUIRE(p.d == 1);
    REQUIRE_FALSE(p.fixed);
    REQUIRE(DegreesOfFreedom<ParticleT>::degreesOfFreedom() == 2);
}

TEST_CASE("CreateTwoParticle") {
    std::mt19937 gen(1);
    auto particles = CreateTwoParticle<ParticleT>::createTwoParticle(gen, areaL);

    REQUIRE(particles.size() == 2);
    REQUIRE((particles[1].r[0] - particles[0].r[0]) == Approx(1.123));

    auto p1 = CreateRandomParticle<ParticleT>::createRandomParticle(gen);
    STATIC_REQUIRE(std::is_same_v<decltype(p1), ParticleT>);
}

TEST_CASE("Getter/Setter Pos/Vel") {
    Vec<Real> pos({3.0, 4.0});
    Vec<Real> vel({3.0, -4.0});
    ParticleT p(0, pos, vel);

    REQUIRE(getGeneralizedPositions(p) == pos); 
    
    Vec<Real> new_pos({9.0, 11.33});
    setGeneralizedPositions(p, new_pos);
    REQUIRE(getGeneralizedPositions(p) == new_pos);

    REQUIRE(getGeneralizedVelocities(p) == vel);
    Vec<Real> new_vel ({-3.2345, 9.6299});
    setGeneralizedVelocities(p, new_vel);
    REQUIRE(getGeneralizedVelocities(p) == new_vel);
}

TEST_CASE("calForceTwo - no PBC adjustment") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}};
    p2.r = {{1.0, 0.0}};

    auto f = calForceTwo(p1, p2, areaL, dummyForce);

    REQUIRE(f[0] == Approx(1.0));
    REQUIRE(f[1] == Approx(0.0));
}

TEST_CASE("calForceTwo - PBC positive correction") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}};
    p2.r = {{6.0, 0.0}};  // dr = (6,0), box=10 => over dr > boxPBC => dr(-4, 0)

    auto f = calForceTwo(p1, p2, areaL, dummyForce);

    REQUIRE(f[0] == Approx(4.0));
    REQUIRE(f[1] == Approx(0.0));
}

TEST_CASE("calForceTwo - PBC negative correction") {
    ParticleT p1, p2;
    p1.r = {{9.0, 0.0}};
    p2.r = {{1.0, 0.0}};  // dr = (-8,0), box=10 => over dr < boxPBC => dr (2, 0)

    auto f = calForceTwo(p1, p2, areaL, dummyForce);

    REQUIRE(f[0] == Approx(2.0));
    REQUIRE(f[1] == Approx(0.0));
}

TEST_CASE("calForceTwo - overlapping particles (r≈0)") {
    ParticleT p1, p2;
    p1.r = {{2.0, 3.0}};
    p2.r = {{2.0, 3.0}}; // the same

    auto f = calForceTwo(p1, p2, areaL, dummyForce);

    REQUIRE(f[0] == Approx(0.0));
    REQUIRE(f[1] == Approx(0.0));
}

TEST_CASE("calTotalForce - single other particle") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}};
    p2.r = {{1.0, 0.0}};

    std::vector<ParticleT> particles{p1, p2};
    auto f = calTotalForce(p1, particles, areaL, dummyForce);

    REQUIRE(f[0] == Approx(1.0));
    REQUIRE(f[1] == Approx(0.0));
}

TEST_CASE("calTotalForce - two other particles, sum up") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}};
    p2.r = {{1.0, 0.0}}; // dr = (1,0)
    p3.r = {{0.0, 2.0}}; // dr = (0,2)

    std::vector<ParticleT> particles{p1, p2, p3};
    auto f = calTotalForce(p1, particles, areaL, dummyForce);

    REQUIRE(f[0] == Approx(1.0));
    REQUIRE(f[1] == Approx(2.0));
}

TEST_CASE("calTotalForce - PBC adjustment applied") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}};
    p2.r = {{6.0, 0.0}}; // dr = 6, box=10 => dr > boxPBC/2 => dr(-4, 0)

    std::vector<ParticleT> particles{p1, p2};
    auto f = calTotalForce(p1, particles, areaL, dummyForce);

    REQUIRE(f[0] == Approx(4.0));
    REQUIRE(f[1] == Approx(0.0));
}

TEST_CASE("calTotalForce - only itself") {
    ParticleT p1;
    p1.r = {{1.0, 1.0}};

    std::vector<ParticleT> particles{p1};
    auto f = calTotalForce(p1, particles, areaL, dummyForce);

    REQUIRE(f[0] == Approx(0.0));
    REQUIRE(f[1] == Approx(0.0));
}

TEST_CASE("calAcceleration - single pair with mass=1") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{1.0, 0.0}}; p2.species = 0;

    std::vector<ParticleT> particles{p1, p2};

    auto a = calAcceleration(p1, particles, species, areaL, dummyForce);

    REQUIRE(a[0] == Approx(1.0 / mass)); // force (1,0) / mass 1
    REQUIRE(a[1] == Approx(0.0 / mass));
}

TEST_CASE("calAcceleration - two neighbors with mass=2") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{2.0, 0.0}}; p2.species = 0; // dr=(2,0)
    p3.r = {{0.0, 4.0}}; p3.species = 0; // dr=(0,4)

    std::vector<ParticleT> particles{p1, p2, p3};
    std::vector<Species<Real>> species{{2.0, 0.0, 0.0}}; // mass=2

    auto a = calAcceleration(p1, particles, species, areaL, dummyForce);

    // totalForce = (2,4), mass=2 → acceleration=(1,2)
    REQUIRE(a[0] == Approx(1.0));
    REQUIRE(a[1] == Approx(2.0));
}

TEST_CASE("calAcceleration - different species mass") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}}; p1.species = 1; // species=1
    p2.r = {{3.0, 0.0}}; p2.species = 0;

    std::vector<ParticleT> particles{p1, p2};
    std::vector<Species<Real>> species{
        {1.0, 0.0, 0.0}, // species 0
        {2.0, 0.0, 0.0}  // species 1
    };

    auto a = calAcceleration(p1, particles, species, areaL, dummyForce);

    // force=(3,0), mass=2 → acceleration=(1.5,0)
    REQUIRE(a[0] == Approx(1.5));
    REQUIRE(a[1] == Approx(0.0));
}

TEST_CASE("calAcceleration - alone particle") {
    ParticleT p1;
    p1.r = {{1.0, 1.0}}; p1.species = 0;

    std::vector<ParticleT> particles{p1};

    auto a = calAcceleration(p1, particles, species, areaL, dummyForce);

    REQUIRE(a[0] == Approx(0.0));
    REQUIRE(a[1] == Approx(0.0));
}

TEST_CASE("calAcceleration - PBC adjustment works") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{6.0, 0.0}}; p2.species = 0; // dr=6, box=10 → adjust to -4

    std::vector<ParticleT> particles{p1, p2};

    auto a = calAcceleration(p1, particles, species, areaL, dummyForce);

    // adjusted dr=(-4,0), force=(-4,0), mass=1 → acc=(-4,0)
    REQUIRE(a[0] == Approx(4.0));
    REQUIRE(a[1] == Approx(0.0));
}

TEST_CASE("calAllAccelerations - single particle") {
    ParticleT p1;
    p1.r = {{0.0, 0.0}};
    p1.species = 0;

    std::vector<ParticleT> particles{p1};

    auto accs = calAllAccelerations(particles, species, areaL, dummyForce);

    REQUIRE(accs.size() == 1);
    REQUIRE(accs[0][0] == Approx(0.0));
    REQUIRE(accs[0][1] == Approx(0.0));
}

TEST_CASE("calAllAccelerations - two symmetric particles") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{2.0, 0.0}}; p2.species = 0;

    std::vector<ParticleT> particles{p1, p2};

    auto accs = calAllAccelerations(particles, species, areaL, dummyForce);

    REQUIRE(accs.size() == 2);
    // For p1: dr from p2=(2,0) → acc=(2,0)
    REQUIRE(accs[0][0] == Approx(2.0));
    REQUIRE(accs[0][1] == Approx(0.0));
    // For p2: dr from p1=(-2,0) → acc=(-2,0)
    REQUIRE(accs[1][0] == Approx(2.0));
    REQUIRE(accs[1][1] == Approx(0.0));
}

TEST_CASE("calAllAccelerations - three particles") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{2.0, 0.0}}; p2.species = 0;
    p3.r = {{0.0, 2.0}}; p3.species = 0;

    std::vector<ParticleT> particles{p1, p2, p3};
    auto accs = calAllAccelerations(particles, species, areaL, dummyForce);

    REQUIRE(accs.size() == 3);
    // For p1: total force=(2,0)+(0,2)=(2,2)
    REQUIRE(accs[0][0] == Approx(2.0));
    REQUIRE(accs[0][1] == Approx(2.0));
    // For p2: total force=(-2,0)+(-2,2)=(-4,2)
    REQUIRE(accs[1][0] == Approx(4.0));
    REQUIRE(accs[1][1] == Approx(2.0));
    // For p3: total force=(0,-2)+(2,-2)=(2,-4)
    REQUIRE(accs[2][0] == Approx(2.0));
    REQUIRE(accs[2][1] == Approx(4.0));
}

TEST_CASE("calAllAccelerations - PBC correction") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{6.0, 0.0}}; p2.species = 0; // dr=6, box=10 → adjust -4

    std::vector<ParticleT> particles{p1, p2};
    std::vector<Species<double>> species{{2.0,0.0,0.0}}; // mass=2

    auto accs = calAllAccelerations(particles, species, areaL, dummyForce);

    REQUIRE(accs.size() == 2);
    // For p1: dr=( -4,0 ), force=(-4,0), mass=2 → acc=(-2,0)
    REQUIRE(accs[0][0] == Approx(2.0));
    REQUIRE(accs[0][1] == Approx(0.0));
    // For p2: dr=( +4,0 ), force=(4,0), mass=2 → acc=(2,0)
    REQUIRE(accs[1][0] == Approx(2.0));
    REQUIRE(accs[1][1] == Approx(0.0));
}

TEST_CASE("implementPBC - inside box remains unchanged") {
    ParticleT p;
    p.r = {{3.0, 7.0}};
    implementPBC(p, areaL);

    REQUIRE(p.r[0] == Approx(3.0));
    REQUIRE(p.r[1] == Approx(7.0));
}

TEST_CASE("implementPBC - values greater than box folded") {
    ParticleT p;
    p.r = {{12.0, 25.0}}; // box=10 → expect (2,5)
    implementPBC(p, areaL);

    REQUIRE(p.r[0] == Approx(2.0));
    REQUIRE(p.r[1] == Approx(5.0));
}

TEST_CASE("implementPBC - negative values corrected") {
    ParticleT p;
    p.r = {{-3.0, -15.0}}; // box=10 → expect (7,5)
    implementPBC(p, areaL);

    REQUIRE(p.r[0] == Approx(7.0));
    REQUIRE(p.r[1] == Approx(5.0));
}

TEST_CASE("implementPBC - multiple box lengths positive") {
    ParticleT p;
    p.r = {{33.0, 47.0}}; // box=10 → expect (3,7)
    implementPBC(p, areaL);

    REQUIRE(p.r[0] == Approx(3.0));
    REQUIRE(p.r[1] == Approx(7.0));
}

TEST_CASE("implementPBC - multiple box lengths negative") {
    ParticleT p;
    p.r = {{-27.0, -42.0}}; // box=10 → expect (3,8)
    implementPBC(p, areaL);

    REQUIRE(p.r[0] == Approx(3.0));
    REQUIRE(p.r[1] == Approx(8.0));
}