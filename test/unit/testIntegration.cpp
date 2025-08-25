#define CATCH_CONFIG_MAIN
#include "./include/catch_amalgamated.hpp"
#include "./include/DummyForce.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/md/Integration.h"

using Catch::Approx;
using Real = double;
using ParticleT = ParticleDot<Real>;
auto dummyForce = DummyForceFunc2D<Real>();

const Real areaL = 10;
const Real mass = 1;
const Real momentI = 1;
const Real radius = 0.5;
std::vector<Species<Real>> species{{mass, momentI, radius}}; // mass=1


// TEST_CASE("EulerSymplecticStep - single free particle no force") {
//     ParticleT p;
//     p.r = {{1.0, 1.0}};
//     p.v = {{0.5, -0.5}};
//     p.species = 0;
//     Real dt = 1.0;

//     std::vector<ParticleT> particles{p};

//     EulerSymplecticStep(particles, species, dt, areaL, dummyForce);

//     // v 不變
//     REQUIRE(particles[0].v[0] == Approx(0.5));
//     REQUIRE(particles[0].v[1] == Approx(-0.5));
//     // r = r + v*dt
//     REQUIRE(particles[0].r[0] == Approx(1.5));
//     REQUIRE(particles[0].r[1] == Approx(0.5));
// }

// TEST_CASE("EulerSymplecticStep - single particle with force") {
//     ParticleT p1, p2;
//     p1.r = {{0.0, 0.0}}; p1.v = {{0.0,0.0}}; p1.species=0;
//     p2.r = {{1.0, 0.0}}; p2.v = {{0.0,0.0}}; p2.species=0;
//     Real dt = 1.0;

//     std::vector<ParticleT> particles{p1, p2};

//     EulerSymplecticStep(particles, species, dt, areaL, dummyForce);

//     // p1: dr from p2=(1,0), so a=(1,0), v=(1,0), r=r+v_old*dt=(0,0)
//     REQUIRE(particles[0].v[0] == Approx(1.0));
//     REQUIRE(particles[0].r[0] == Approx(0.0));

//     // p2: dr from p1=(-1,0), so a=(-1,0), v=(-1,0), r=r+v_old*dt=(1,0)
//     REQUIRE(particles[1].v[0] == Approx(-1.0));
//     REQUIRE(particles[1].r[0] == Approx(1.0));
// }

// TEST_CASE("EulerSymplecticStep - fixed particle remains unchanged") {
//     ParticleT p;
//     p.r = {{5.0, 5.0}};
//     p.v = {{1.0, 1.0}};
//     p.fixed = true;
//     p.species = 0;

//     std::vector<ParticleT> particles{p};
//     std::vector<Species<double>> species{{1.0,0.0,0.0}};

//     EulerSymplecticStep(particles, species, 1.0, 10.0, dummyForce);

//     // r,v unchanged
//     REQUIRE(particles[0].r[0] == Approx(5.0));
//     REQUIRE(particles[0].r[1] == Approx(5.0));
//     REQUIRE(particles[0].v[0] == Approx(1.0));
//     REQUIRE(particles[0].v[1] == Approx(1.0));
// }

// TEST_CASE("EulerSymplecticStep - PBC correction applied") {
//     ParticleT p1, p2;
//     p1.r = {{9.5, 0.0}}; p1.v = {{1.0,0.0}}; p1.species=0;
//     p2.r = {{0.0, 0.0}}; p2.v = {{0.0,0.0}}; p2.species=0;

//     std::vector<ParticleT> particles{p1, p2};
//     std::vector<Species<double>> species{{1.0,0.0,0.0}};

//     EulerSymplecticStep(particles, species, 1.0, 10.0, [](auto&,auto&,auto&,auto){return Vec<double,2>{{0,0}};});

//     // 原本 r=9.5, v=1, dt=1 → r+v*dt=10.5 → fold to 0.5
//     REQUIRE(particles[0].r[0] == Approx(0.5));
// }