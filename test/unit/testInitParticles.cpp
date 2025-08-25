#include <vector>
#include <ostream>
#include <iostream>
#include <string>
#include <random>
#include <cassert>

#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/CircleDistribution.h"

using Real = double;
using ParticleT = ParticleDot<Real>;

const Real mass = 1;
const Real momentI = 1;
const Real radius = 0.5;
unsigned int particlesNum = 49;
const Real areaL = 20;
const Real boxPBC = 20;
const Real fixRadius = 20;
int speciesInd = 0;
Species<Real> species1 {mass, momentI, radius};
Species<Real> species2 {2.0f * mass, momentI, 0.5f * radius};
std::vector<Species<Real>> allSpecies {species1, species2};
const std::string particlesInit = "RANDOM";
const std::string particlesType = "RRUU";


std::vector<std::vector<Real>> default49_particles {
    {1.94219, 3.24856}, {2.33187, 5.25462}, {2.21991, 7.53682}, {2.55162, 10.6209},
    {3.15869, 12.9176}, {1.8892, 15.0272}, {2.99441, 17.1596}, {5.1396, 2.75748},
    {4.68444, 4.46318}, {4.30126, 7.68604}, {4.9192, 9.58319}, {4.39426, 13.1051},
    {4.37543, 15.6253}, {5.12658, 17.5989}, {7.14147, 2.09652}, {7.48959, 5.0609},
    {6.95995, 7.93861}, {7.56118, 9.92197}, {7.55764, 12.7284}, {7.70675, 14.4395},
    {7.28091, 17.8948}, {10.572, 1.7675}, {10.4304, 4.3461}, {9.8197, 7.89438},
    {10.4091, 9.47939}, {10.3136, 12.456}, {9.6535, 15.4975}, {9.95871, 17.8641},
    {12.0955, 2.51237}, {12.5283, 4.50842}, {12.5147, 6.9032}, {13.2015, 9.66748},
    {12.6051, 12.3769}, {12.1846, 15.3261}, {12.8778, 17.3917}, {14.48, 3.15579},
    {15.1745, 4.304}, {14.6915, 7.29122}, {14.4782, 9.49243}, {14.325, 12.8614},
    {15.7226, 15.5895}, {14.3095, 17.9708}, {17.5477, 1.914}, {17.0735, 5.62467},
    {17.8964, 7.0691}, {16.9933, 9.7847}, {18.1465, 12.6149}, {17.2368, 15.342},
    {16.9981, 17.7806}
};

template<typename TParticle, typename T = typename TParticle::value_type>
std::ostream& operator<<(std::ostream& out, const TParticle& p) {
    auto pos = getGeneralizedPositions(p);
    out << "( " << pos[0] << ", " << pos[1] << ")";
    return out;
}


std::mt19937 GetRandomGenerator(unsigned int seed=1) {
    if (seed > 0) {
        return std::mt19937(seed);
    } else {
        std::random_device rd;
        return std::mt19937(rd());
    }
}

const Real EPSILON = 1e-3;
void Test_RANDOM() {
    auto gen = GetRandomGenerator();
    auto particles = initialParticles<ParticleT>(particlesNum, allSpecies, speciesInd, areaL, gen, particlesInit, particlesType);
    assert(particles.size() == particlesNum);
    
    assert(default49_particles.size() == particles.size());
    for (int i = 0; i < particles.size(); i++) {
        auto pos1 = getGeneralizedPositions(particles[i]);
        auto pos2 = default49_particles[i];
        // std::cerr << "pos1[0] = " << pos1[0] << ", pos2[0] = " << pos2[0] << std::endl;
        // std::cerr << "pos1[1] = " << pos1[1] << ", pos2[1] = " << pos2[1] << std::endl;
        assert(std::abs(pos1[0] - pos2[0]) < EPSILON);
        assert(std::abs(pos1[1] - pos2[1]) < EPSILON);
    }
        
}

int main() {
    Test_RANDOM();

    std::cout << "All tests passed!\n";
    return 0;
}
