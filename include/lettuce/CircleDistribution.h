#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <sstream>
#include <stdexcept>
#include <adios2.h>

#include "Vec.h"
#include "Circle.h"
#include "lettuce/ParticleDot.h"
#include "lettuce/ParticleOriented.h"
#include "lettuce/CirclesIntersectionFuncs.h"

// ================================== Circle-Particle conversion ==================================
template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<TParticle> circleToParticle(const std::vector<Circle<T>> &circles) {
    std::vector<TParticle> particles;
    
    for (auto &c: circles){
        TParticle p;
        p.r = c.c;
        particles.push_back(p);
    }
    return particles;
}

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<Circle<T>> particleToCircle(const std::vector<TParticle> &particles, const std::vector<Species<T>>& allSpecies, const int& speciesNum) {
    std::vector<Circle<T>> circles;
    Circle<T> c;

    for (auto &p: particles){
        c.c = p.r;
        c.r = allSpecies[speciesNum].radius;
        circles.push_back(c);
    }
    return circles;
}

// ================================== old Random Particle ==================================

template<typename T>
std::vector<Circle<T>> distRandomCircles(const int& circlesNum, const T& L, const T& radius, std::mt19937& gen) {
    std::vector<Circle<T>> circles;
    std::uniform_real_distribution<> dis(0, L);
    
    while (circles.size() < circlesNum) {
        auto newCircle = placeRandomCircle(circles, dis, radius, gen);
        if (!std::isnan(newCircle.c[0])) {
            circles.push_back(newCircle);
        }
    }
    return circles;
}

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<TParticle> distRandomParticles(const int& particlesNum, const T& L, const T& radius, std::mt19937& gen) {
    std::vector<TParticle> randomParticles(particlesNum);
    std::vector<Circle<T>> randomCircles = distRandomCircles(particlesNum, L, radius, gen);
    for (int i = 0; i < particlesNum; ++i) {
        randomParticles[i].r = randomCircles[i].c;
    }
    return randomParticles;
}

template <typename T>
bool has_overlap(const Circle<T> &c1, const Circle<T> &c2) {
    T distance2 = (c1.c - c2.c).abs2();
    const auto radiuses = c1.r + c2.r;
    return distance2 < radiuses * radiuses;
}

template<typename T>
bool has_overlap(const Circle<T>& newCircle, const std::vector<Circle<T>>& circles) {
    for (const auto& circle : circles) {
        if (has_overlap(newCircle, circle)) {
            return true;
        }
    }
    return false;
}

template<typename T>
Circle<T> placeRandomCircle(const std::vector<Circle<T>>& circles, std::uniform_real_distribution<> dis, const T& radius, std::mt19937& gen) {
    Circle<T> newCircle;
    newCircle.c[0] = dis(gen);
    newCircle.c[1] = dis(gen);
    newCircle.r = radius;

    if (has_overlap(newCircle, circles)) {
        newCircle.c[0] = NAN;  // Indicate overlap
    }
    return newCircle;
}

// ================================== Manual Random Particle ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<TParticle> manualRandomParticles(const int& particlesNum, T L, const T& radius, std::mt19937& gen) {
    constexpr double maxPackingDensity = 0.7854;

    T requiredAreaL = std::sqrt(particlesNum * M_PI * radius * radius / maxPackingDensity);

    if (L < requiredAreaL) {
        std::cout << "The specified area (" << L << ") is too small for " << particlesNum 
                  << " particles with radius " << radius << ".\n";
        std::cout << "Minimum required areaL: " << requiredAreaL << std::endl;
        throw std::runtime_error("Insufficient area for particle placement.");
    }

    int topSquareRoot = static_cast<int>(std::ceil(std::sqrt(particlesNum)));
    std::vector<TParticle> particles;
    particles.reserve(particlesNum);

    T distance = L / (topSquareRoot + 1);
    const T dr = distance / 2 - radius;
    std::uniform_real_distribution<T> randomDisplacement(-dr, dr);

    for (int i = 0; i < topSquareRoot; ++i) {
        for (int j = 0; j < topSquareRoot; ++j) {
            if (particles.size() < particlesNum) {
                auto newParticle = createRandomParticle<TParticle>(gen);
                newParticle.r[0] = (i + 1) * distance + randomDisplacement(gen);
                newParticle.r[1] = (j + 1) * distance + randomDisplacement(gen);

                bool overlap = false;
                for (const auto& existingParticle : particles) {
                    T dx = newParticle.r[0] - existingParticle.r[0];
                    T dy = newParticle.r[1] - existingParticle.r[1];
                    if ((dx * dx + dy * dy) < (4 * radius * radius)) {
                        overlap = true;
                        break;
                    }
                }
                if (!overlap) {
                    particles.push_back(newParticle);
                } else {
                    --j;
                }
            }
        }
    }

    return particles;
}

// ================================== DLA shoot Particle ==================================

template<typename T>
Circle<T> startCircleRandom(const T& radius, const T& areaRadius, std::mt19937 &gen) {   
    Circle<T> newCircle;
    newCircle.r = radius;

    std::uniform_real_distribution<> randomAngle(0, 2 * M_PI);
    T angle = randomAngle(gen);

    newCircle.c[0] = areaRadius * cos(angle);
    newCircle.c[1] = areaRadius * sin(angle);

    return newCircle;
}

template<typename T>
Vec<T> shootToCenter(const Circle<T> &startCircle, const T &areaWidth) {
    return (-startCircle.c + areaWidth / static_cast<T>(2));
}

template<typename T>
Circle<T> findStopPointAll (const Circle<T> &startCircle, Vec<T> &direction, const std::vector<Circle<T>> &circles){
    std::vector<std::pair<Circle<T>, T>> pairs;
    for (const auto &c: circles){
        pairs.push_back(findStopPoint(startCircle, direction, c));
    }

    if (pairs.empty()) {
        return startCircle;
    }

    auto it = std::min_element(pairs.begin(), pairs.end(),
        [](const std::pair<Circle<T>, T>& a, const std::pair<Circle<T>, T>& b) {
            if (!std::isnan(a.second) && !std::isnan(b.second)) {
                return a.second < b.second;
            }
            else if (std::isnan(a.second)) {
                return false;
            }
            else {
                return true;
            }
        });

    return it->first;
}

template<typename T>
bool isWithinBounds(const Circle<T>& circle, const T& L) {
    return (circle.c[0] - circle.r >= 0 && circle.c[0] + circle.r <= L &&
            circle.c[1] - circle.r >= 0 && circle.c[1] + circle.r <= L);
}

template<typename T>
std::vector<Circle<T>> distCirclesDLA(const int& shootNum, const T& L, const T& radius, std::mt19937& gen) {
    const double packingDensity = 0.70;
    T area = L * L;
    T circleArea = M_PI * radius * radius;
    int maxCircles = static_cast<int>(packingDensity * area / circleArea);

    Circle<T> target = {{{L / 2, L / 2}}, radius};
    std::vector<Circle<T>> finalCircles;
    finalCircles.push_back(target);

    while (finalCircles.size() < shootNum -1 && finalCircles.size() < maxCircles) {
        Circle<T> newCircle = startCircleRandom(radius, L * 10, gen);
        Vec<T> direction = shootToCenter(newCircle, L);
        Circle<T> endPoint = findStopPointAll(newCircle, direction, finalCircles);
        if (!std::isnan(endPoint.c[0]) && isWithinBounds(endPoint, L)) {
            finalCircles.push_back(endPoint);
        }
    }

    if (finalCircles.size() < shootNum - 2) {
        std::cout << "With the density of: " << packingDensity * 100 << " %, only " << finalCircles.size() << " circles could be placed." << std::endl;
    }

    return finalCircles;
}

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<TParticle> distParticleDLA(const int& shootNum, const T& L, const T& radius, std::mt19937& gen) {
    std::vector<Circle<T>> finalCircles = distCirclesDLA(shootNum, L, radius, gen);
    std::vector<TParticle> finalParticles(finalCircles.size());

    for (size_t i = 0; i < finalCircles.size(); ++i) {
        finalParticles[i] = createRandomParticle<TParticle>(gen);
        finalParticles[i].r = finalCircles[i].c;
    }
    return finalParticles;
}

// ================================== Initial all Particles ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<TParticle> initialParticles(const unsigned int& particlesNum, 
                                        const std::vector<Species<T>>& allSpecies, 
                                        const int& speciesNum, 
                                        const T& L, std::mt19937& gen, 
                                        const std::string& configuration) {
    std::vector<TParticle> particles;

    if (configuration == "RANDOM") {
        particles = manualRandomParticles<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } else if (configuration == "RANDOM2") {
        particles = distRandomParticles<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } else if (configuration == "DLA") {
        particles = distParticleDLA<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } else if (configuration == "TWO") {
        particles = createTwoParticle<TParticle, T>(gen, L);
    }
    else if (configuration.ends_with(".bp")) {

        adios2::ADIOS adios;
        adios2::IO io = adios.DeclareIO("ReadConfig");

        adios2::Engine engine = io.Open(configuration, adios2::Mode::ReadRandomAccess);

        adios2::Variable<T> varPositions = io.InquireVariable<T>("positions");
        adios2::Variable<T> varVelocities = io.InquireVariable<T>("velocities");
        if (!varPositions || !varVelocities) {
            throw std::runtime_error("Missing required variables in file: " + configuration);
        }

        size_t totalSteps = varPositions.Steps();
        size_t lastStep = totalSteps - 1;
        // std::cout << "Reading step " << lastStep + 1 << " (final step) of " << totalSteps << " steps." << std::endl;
        
        varPositions.SetStepSelection({lastStep, 1});
        varVelocities.SetStepSelection({lastStep, 1});

        std::vector<size_t> shape = varPositions.Shape();
        const size_t num_particles = shape[0];
        const size_t dimensions = shape[1];

        std::vector<T> positions(num_particles * dimensions);
        std::vector<T> velocities(num_particles * dimensions);

        engine.Get(varPositions, positions.data(), adios2::Mode::Sync);
        engine.Get(varVelocities, velocities.data(), adios2::Mode::Sync);
        engine.Close();

        for (size_t i = 0; i < num_particles; ++i) {
            TParticle p;
            constexpr int D = degreesOfFreedom<TParticle>();
            
            Vec<T, D> q, qDot;
            for (int a = 0; a < D; ++a) {
                q[a] = positions[i * dimensions + a];
                qDot[a] = velocities[i * dimensions + a];
            }

            if (q[0] <= L && q[0] >= 0. && q[1] <= L && q[1] >= 0.) {
                setGeneralizedPositions(p, q);
                setGeneralizedVelocities(p, qDot);
                particles.push_back(p);
            }
        }
    }

    // Set species for all particles
    for (auto& p : particles) {
        p.species = speciesNum;
    }

    return particles;
}

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<ParticleOriented<T>> initialParticlesOriented(const unsigned int& particlesNum, 
                                                          const std::vector<Species<T>>& allSpecies, 
                                                          const int& speciesNum, 
                                                          const T& L, std::mt19937& gen, 
                                                          const std::string& configuration) {
    std::vector<ParticleOriented<T>> particlesOriented;
    std::vector<TParticle> particlesDot = initialParticles<TParticle>(particlesNum, allSpecies, speciesNum, L, gen, configuration);

    particlesOriented.resize(particlesDot.size());
    
    for (size_t i = 0; i < particlesDot.size(); i++) {
        particlesOriented[i].r = particlesDot[i].r;
        particlesOriented[i].v = particlesDot[i].v;
        particlesOriented[i].species = particlesDot[i].species;
        particlesOriented[i].phi = 0;
    }
    
    return particlesOriented;
}

// ================================== Deposite Particles RANDOM ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
void addParticle(std::vector<TParticle>& particles,  
                 const std::vector<Species<T>>& allSpecies, 
                 const int& speciesNum, 
                 const T& L, 
                 const std::string& depositeMethod, 
                 std::mt19937& gen) 
{
    int particlesNum = 1;
    TParticle newParticle;

    do {            
        if (depositeMethod == "RANDOM") {
            auto tempParticles = manualRandomParticles<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
            if (!tempParticles.empty()) newParticle = tempParticles[0];
        } else if (depositeMethod == "RANDOM2") {
            auto tempParticles = distRandomParticles<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
            if (!tempParticles.empty()) newParticle = tempParticles[0];
        } else if (depositeMethod == "DLA") {
            auto tempParticles = distParticleDLA<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
            if (!tempParticles.empty()) newParticle = tempParticles[0];
        } else {
            throw std::runtime_error("Unknown deposition method: " + depositeMethod);
        }
    } while (hasOverlap(newParticle, particles, allSpecies, speciesNum));

    particles.push_back(newParticle);
}

template<typename TParticle, typename T = typename TParticle::value_type>
bool hasOverlap(const TParticle &p1, const TParticle &p2, const std::vector<Species<T>>& allSpecies, const int& speciesNum) {
    T distance2 = (p1.r - p2.r).abs2();
    const auto radiuses = allSpecies[speciesNum].radius * 2;
    return distance2 < radiuses * radiuses;
}

template<typename TParticle, typename T = typename TParticle::value_type>
bool hasOverlap(const TParticle &newParticle, const std::vector<TParticle> particles, const std::vector<Species<T>>& allSpecies, const int& speciesNum) {
    for (const auto& particle : particles) {
        if (hasOverlap(newParticle, particle, allSpecies, speciesNum)) {
            return true;
        }
    }
    return false;
}

// ================================== Deposite Particles DLA ==================================

template<typename TParticle, typename T = typename TParticle::value_type>
auto depositeDLAinfo(const std::vector<TParticle>& particles,  
                         const std::vector<Species<T>>& allSpecies, 
                         const int& speciesNum, 
                         const T& L,
                         std::mt19937& gen) 
{
    std::vector<Circle<T>> circles = particleToCircle<TParticle, T>(particles, allSpecies, speciesNum);

    T radius = allSpecies[speciesNum].radius;
    Circle<T> endPoint;

    int attempt = 0;
    int maxAttempt = 10000000;
    
    // std::cout << ">>>  shooting new particle:" << attempt << std::endl;
    while (!isWithinBounds(endPoint, L) && attempt<maxAttempt) {
        Circle<T> newCircle = startCircleRandom(radius, L * 10, gen);
        Vec<T> direction = shootToCenter(newCircle, L);
        endPoint = findStopPointAll(newCircle, direction, circles);
        attempt ++;
    }
    // std::cout << ">> attempt numbers for successful DLA: " << attempt << std::endl;

    if (std::isnan(endPoint.c[0])){
        std::cout << "could not add new particle!!" << std::endl;
    }

    return endPoint.c;
}

// ================================== Useless ones (clean later) ==================================

template<typename T>
T packingDensity(const int& number, const T& radius, const T& L) {
    return (number * M_PI * radius * radius) / (L * L);
}

template<typename T>
void writeCircles(T begin, T end, const std::string& fname) {
    std::ofstream output_file(fname);
    for (auto c = begin; c != end; ++c) {
        output_file << c->c << ", " << c->r << "\n";
    }
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeParticle(const std::vector<TParticle>& particles, T radius, const std::string& fname) {
    std::ofstream output_file(fname);
    for (const auto& p : particles) {
        output_file << p.r[0] << ", " << p.r[1] << ", " << radius << "\n";
    }
}
