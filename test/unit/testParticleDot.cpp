#define CATCH_CONFIG_MAIN
#include <random>
#include "./include/catch_amalgamated.hpp"
#include "./include/DummyForce.h"
#include "lettuce/core/Vec.h"
#include "lettuce/core/Species.h"
#include "lettuce/core/ParticleDot.h"

using Catch::Approx;
using Real = double;
using ParticleT = ParticleDot<Real>;
// auto dummyForce = DummyForceFunc2D<Real>();
const Real LJepsilon = 1;
const Real LJsigma = 1;
auto ForceFunc = LennardJonesForce2D<Real>(LJepsilon, LJsigma);
const Real areaL = 10;
const Real mass = 1;
const Real momentI = 1;
const Real radius = 0.5;
std::vector<Species<Real>> species{{mass, momentI, radius}}; // mass=1


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

    auto f_p1 = calForceTwo(p1, p2, areaL, ForceFunc);
    auto f_p2 = calForceTwo(p2, p1, areaL, ForceFunc);

    REQUIRE(f_p1.abs() == f_p2.abs());
}

TEST_CASE("calForceTwo - PBC positive correction with norm check") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}};
    p2.r = {{6.0, 0.0}};  // dr = (6,0), after PBC => (-4,0)

    auto f_p1 = calForceTwo(p1, p2, areaL, ForceFunc);
    auto f_p2 = calForceTwo(p2, p1, areaL, ForceFunc);

    // Magnitudes should match (||F12|| = ||F21||)
    REQUIRE(f_p1.abs() == Approx(f_p2.abs()));
}


TEST_CASE("calForceTwo - PBC negative correction") {
    ParticleT p1, p2;
    p1.r = {{9.0, 0.0}};
    p2.r = {{1.0, 0.0}};  // dr = (-8,0), after PBC => (2,0)

    auto f_p1 = calForceTwo(p1, p2, areaL, ForceFunc); // F on p1 due to p2
    auto f_p2 = calForceTwo(p2, p1, areaL, ForceFunc); // F on p2 due to p1

    REQUIRE(f_p1.abs() == f_p2.abs());
}

TEST_CASE("calForceTwo - overlapping particles (r≈0) with LJForceFunc") {
    ParticleT p1, p2;
    p1.r = {{2.0, 3.0}};
    p2.r = {{2.0, 3.0}}; // same position → r=0

    // Force on p1 due to p2
    auto f_p1 = calForceTwo(p1, p2, areaL, ForceFunc);
    // Force on p2 due to p1
    auto f_p2 = calForceTwo(p2, p1, areaL, ForceFunc);

    // Expect zero vector (because r≈0 handled explicitly)
    REQUIRE(f_p1.abs() == f_p2.abs());
    REQUIRE(f_p1.abs() == Approx(0.0));
    // Newton's 3rd law still holds (0 == -0)
    REQUIRE(f_p1[0] == Approx(-f_p2[0]));
    REQUIRE(f_p1[1] == Approx(-f_p2[1]));
}


TEST_CASE("calTotalForce - three particles in x-axis") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}};
    p2.r = {{1.0, 0.0}};
    p3.r = {{2.0, 0.0}};

    std::vector<ParticleT> particles{p1, p2, p3};

    auto f_p1 = calTotalForce(p1, particles, areaL, ForceFunc);
    auto f_p2 = calTotalForce(p2, particles, areaL, ForceFunc);
    auto f_p3 = calTotalForce(p3, particles, areaL, ForceFunc);

    auto f12 = calForceTwo(p1, p2, areaL, ForceFunc);
    auto f21 = calForceTwo(p2, p1, areaL, ForceFunc);
    REQUIRE(f12.abs() == Approx(f21.abs()));

    auto f23 = calForceTwo(p2, p3, areaL, ForceFunc);
    auto f32 = calForceTwo(p3, p2, areaL, ForceFunc);
    REQUIRE(f23.abs() == Approx(f32.abs()));

    auto f13 = calForceTwo(p1, p3, areaL, ForceFunc);
    auto f31 = calForceTwo(p3, p1, areaL, ForceFunc);
    REQUIRE(f13.abs() == Approx(f31.abs()));

    auto f_total = f_p1 + f_p2 + f_p3;
    REQUIRE(f_total.abs() == Approx(0.0));
}

TEST_CASE("calTotalForce - three particles in y-axis") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}};
    p2.r = {{0.0, 6.0}};
    p3.r = {{0.0, 10.0}};

    std::vector<ParticleT> particles{p1, p2, p3};

    auto f_p1 = calTotalForce(p1, particles, areaL, ForceFunc);
    auto f_p2 = calTotalForce(p2, particles, areaL, ForceFunc);
    auto f_p3 = calTotalForce(p3, particles, areaL, ForceFunc);

    auto f12 = calForceTwo(p1, p2, areaL, ForceFunc);
    auto f21 = calForceTwo(p2, p1, areaL, ForceFunc);
    REQUIRE(f12.abs() == Approx(f21.abs()));

    auto f23 = calForceTwo(p2, p3, areaL, ForceFunc);
    auto f32 = calForceTwo(p3, p2, areaL, ForceFunc);
    REQUIRE(f23.abs() == Approx(f32.abs()));

    auto f13 = calForceTwo(p1, p3, areaL, ForceFunc);
    auto f31 = calForceTwo(p3, p1, areaL, ForceFunc);
    REQUIRE(f13.abs() == Approx(f31.abs()));

    auto f_total = f_p1 + f_p2 + f_p3;
    REQUIRE(f_total.abs() == Approx(0.0));
}


TEST_CASE("calTotalForce - two other particles, sum up with norm and 3rd law") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}};
    p2.r = {{1.0, 0.0}}; // dr = (1,0)
    p3.r = {{0.0, 2.0}}; // dr = (0,2)

    std::vector<ParticleT> particles{p1, p2, p3};

    auto f_p1 = calTotalForce(p1, particles, areaL, ForceFunc);
    auto f_p2 = calTotalForce(p2, particles, areaL, ForceFunc);
    auto f_p3 = calTotalForce(p3, particles, areaL, ForceFunc);

    // pairwise forces (Newton's 3rd law)
    auto f12 = calForceTwo(p1, p2, areaL, ForceFunc);
    auto f21 = calForceTwo(p2, p1, areaL, ForceFunc);
    REQUIRE(f12.abs() == Approx(f21.abs()));

    auto f23 = calForceTwo(p2, p3, areaL, ForceFunc);
    auto f32 = calForceTwo(p3, p2, areaL, ForceFunc);
    REQUIRE(f23.abs() == Approx(f32.abs()));

    auto f13 = calForceTwo(p1, p3, areaL, ForceFunc);
    auto f31 = calForceTwo(p3, p1, areaL, ForceFunc);
    REQUIRE(f13.abs() == Approx(f31.abs()));

    auto f_total = f_p1 + f_p2 + f_p3;
    REQUIRE(f_total.abs() == Approx(0.0).margin(1e-12));
}


TEST_CASE("calTotalForce - PBC adjustment applied with ForceFunc") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}};
    p2.r = {{6.0, 0.0}}; // dr=6, box=10 => adjust to (-4,0)

    std::vector<ParticleT> particles{p1, p2};
    auto f_p1 = calTotalForce(p1, particles, areaL, ForceFunc);
    auto f_p2 = calTotalForce(p2, particles, areaL, ForceFunc);

    // pairwise check
    auto f12 = calForceTwo(p1, p2, areaL, ForceFunc);
    auto f21 = calForceTwo(p2, p1, areaL, ForceFunc);

    REQUIRE(f12.abs() == Approx(f21.abs()));

    // total force = 0
    auto f_total = f_p1 + f_p2;
    REQUIRE(f_total.abs() == Approx(0.0).margin(1e-12));
}

TEST_CASE("calAcceleration - single pair with mass=1") {
    ParticleT p1, p2;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{7.0, 0.0}}; p2.species = 0;

    std::vector<ParticleT> particles{p1, p2};

    auto acc_p1 = calAcceleration(p1, particles, species, areaL, ForceFunc);
    auto acc_p2 = calAcceleration(p2, particles, species, areaL, ForceFunc);

    // total force = 0
    auto acc_total = acc_p1 + acc_p2;
    REQUIRE(acc_total.abs() == Approx(0.0).margin(1e-12));
}

TEST_CASE("calAcceleration - two neighbors with mass=2") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{2.0, 0.0}}; p2.species = 0; // dr=(2,0)
    p3.r = {{0.0, 4.0}}; p3.species = 0; // dr=(0,4)

    std::vector<ParticleT> particles{p1, p2, p3};
    std::vector<Species<Real>> species{{2.0, 0.0, 0.0}}; // mass=2

    auto acc_p1 = calAcceleration(p1, particles, species, areaL, ForceFunc);
    auto acc_p2 = calAcceleration(p2, particles, species, areaL, ForceFunc);
    auto acc_p3 = calAcceleration(p3, particles, species, areaL, ForceFunc);

    auto acc_total = acc_p1 + acc_p2 + acc_p3;
    REQUIRE(acc_total.abs() == Approx(0.0).margin(1e-12));
}

TEST_CASE("calAcceleration - alone particle") {
    ParticleT p1;
    p1.r = {{1.0, 1.0}}; p1.species = 0;

    std::vector<ParticleT> particles{p1};

    auto a = calAcceleration(p1, particles, species, areaL, ForceFunc);

    REQUIRE(a[0] == Approx(0.0));
    REQUIRE(a[1] == Approx(0.0));
}

TEST_CASE("calAcceleration - x-axis PBC adjustment works") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{6.0, 0.0}}; p2.species = 0;
    p3.r = {{9.0, 0.0}}; p3.species = 0;

    std::vector<ParticleT> particles{p1, p2, p3};

    auto acc_p1 = calAcceleration(p1, particles, species, areaL, ForceFunc);
    auto acc_p2 = calAcceleration(p2, particles, species, areaL, ForceFunc);
    auto acc_p3 = calAcceleration(p3, particles, species, areaL, ForceFunc);

    auto acc_total = acc_p1 + acc_p2 + acc_p3;
    REQUIRE(acc_total.abs() == Approx(0.0).margin(1e-12));
}

TEST_CASE("calAcceleration - y-axis PBC adjustment works") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{0.0, 6.0}}; p2.species = 0;
    p3.r = {{0.0, 9.0}}; p3.species = 0;

    std::vector<ParticleT> particles{p1, p2, p3};

    auto acc_p1 = calAcceleration(p1, particles, species, areaL, ForceFunc);
    auto acc_p2 = calAcceleration(p2, particles, species, areaL, ForceFunc);
    auto acc_p3 = calAcceleration(p3, particles, species, areaL, ForceFunc);

    auto acc_total = acc_p1 + acc_p2 + acc_p3;
    REQUIRE(acc_total.abs() == Approx(0.0).margin(1e-12));
}


TEST_CASE("calAllAccelerations - 3 particles") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{2.0, 0.0}}; p2.species = 0; // dr=(2,0)
    p3.r = {{0.0, 4.0}}; p3.species = 0; // dr=(0,4)

    std::vector<ParticleT> particles{p1, p2, p3};

    auto cal_all_acc = calAllAccelerations(particles, species, areaL, ForceFunc);

    
    Vec<Real, 2> sum;
    for (const auto & each_acc: cal_all_acc) {
        sum += each_acc;
    }
    REQUIRE(sum.abs() == Approx(0.0).margin(1e-12));
}

TEST_CASE("calAllAccelerations - 3 particles, along with x-axis") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{6.0, 0.0}}; p2.species = 0;
    p3.r = {{9.0, 0.0}}; p3.species = 0;

    std::vector<ParticleT> particles{p1, p2, p3};

    auto cal_all_acc = calAllAccelerations(particles, species, areaL, ForceFunc);

    Vec<Real, 2> sum;
    for (const auto & each_acc: cal_all_acc) {
        sum += each_acc;
    }
    REQUIRE(sum.abs() == Approx(0.0).margin(1e-12));
}

TEST_CASE("calAllAccelerations - 3 particles, along with y-axis") {
    ParticleT p1, p2, p3;
    p1.r = {{0.0, 0.0}}; p1.species = 0;
    p2.r = {{0.0, 6.0}}; p2.species = 0;
    p3.r = {{0.0, 9.0}}; p3.species = 0;

    std::vector<ParticleT> particles{p1, p2, p3};

    auto cal_all_acc = calAllAccelerations(particles, species, areaL, ForceFunc);

    Vec<Real, 2> sum;
    for (const auto & each_acc: cal_all_acc) {
        sum += each_acc;
    }
    REQUIRE(sum.abs() == Approx(0.0).margin(1e-12));
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