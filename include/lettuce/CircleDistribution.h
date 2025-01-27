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
bool isWithinBounds(const Circle<T>& circle, const T& L) {
    return (circle.c[0] - circle.r >= 0 && circle.c[0] + circle.r <= L &&
            circle.c[1] - circle.r >= 0 && circle.c[1] + circle.r <= L);
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

template<typename T>
T packingDensity(const int& number, const T& radius, const T& L) {
    return (number * M_PI * radius * radius) / (L * L);
}

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

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<TParticle> initialParticles(const unsigned int& particlesNum, const std::vector<Species<T>>& allSpecies, int& speciesNum, const T& L, std::mt19937& gen, const std::string& configuration) {
    std::vector<TParticle> particles;

    if (configuration == "RANDOM") {
        particles = manualRandomParticles<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } else if (configuration == "RANDOM2") {
        particles = distRandomParticles<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } else if (configuration == "DLA") {
        particles = distParticleDLA<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } 
    else if (configuration.ends_with(".bp")) {
        particles.reserve(particlesNum);
        std::cout << "Reading configuration from file: " << configuration << std::endl;
        adios2::ADIOS adios;
        adios2::IO io = adios.DeclareIO("ReadConfig");
        adios2::Engine engine = io.Open(configuration, adios2::Mode::Read);

        engine.BeginStep();
        adios2::Variable<T> varPositions = io.InquireVariable<T>("positions");
        adios2::Variable<T> varVelocities = io.InquireVariable<T>("velocities");

        if (!varPositions || !varVelocities) {
            throw std::runtime_error("Missing 'positions' or 'velocities' in ADIOS2 file: " + configuration);
        }

        size_t totalSteps = varPositions.Steps();
        if (totalSteps == 0) {
            throw std::runtime_error("No steps found in ADIOS2 file.");
        }
        const size_t lastStep = totalSteps - 1;

        varPositions.SetStepSelection({lastStep, 1});
        varVelocities.SetStepSelection({lastStep, 1});

        std::vector<size_t> shape = varPositions.Shape();
        if (shape[0] != particlesNum || shape[1] != degreesOfFreedom<TParticle>()) {
            throw std::runtime_error("Shape of positions in ADIOS file does not match expected number of particles or degrees of freedom.");
        }

        std::vector<T> positions(shape[0] * shape[1]);
        std::vector<T> velocities(shape[0] * shape[1]);

        engine.Get(varPositions, positions.data());
        engine.Get(varVelocities, velocities.data());
        engine.PerformGets();
        engine.EndStep();
        engine.Close();

        for (size_t i = 0; i < shape[0]; ++i) {
            TParticle p;
            constexpr int D = degreesOfFreedom<TParticle>();
            Vec<T, D> q, qDot;
            for (int a = 0; a < shape[1]; ++a) {
                q[a] = positions[i * shape[1] + a];
                qDot[a] = velocities[i * shape[1] + a];
            }

            setGeneralizedPositions(p, q);
            setGeneralizedVelocities(p, qDot);
            particles.push_back(p);
        }
    }  
    else { // Old .dat format
        particles.reserve(1024);
        std::ifstream config(configuration);

        if (!config.is_open()) {
            throw std::runtime_error("Could not open configuration file: " + configuration);
        }
        std::string line;
        while (std::getline(config, line)) {
            TParticle p;
            std::istringstream is(line);
            constexpr int D = degreesOfFreedom<TParticle>();
            // Read positions
            Vec<T, D> q;
            for (int a = 0; a < D; ++a) {
                if (!(is >> q[a])) {
                    std::ostringstream os;
                    os << "Invalid particle position line " << particles.size() + 1;
                    throw std::runtime_error(std::move(os).str());
                }
            }
            setGeneralizedPositions(p, q);
            // Filter out particle outside of LxL box
            if (q[0] > L || q[0] < 0. || q[1] > L || q[1] < 0.) {
                continue;
            }
            // Read velocities
            Vec<T, D> qDot;
            for (int a = 0; a < D; ++a) {
                if (!(is >> qDot[a])) {
                    qDot.fill(0);
                    break;
                }
            }
            setGeneralizedVelocities(p, qDot);
            if (!is.eof()) {
                std::ostringstream os;
                os << "Invalid particle velocity line " << particles.size() + 1 << ": unread characters.";
                throw std::runtime_error(std::move(os).str());
            }
            particles.push_back(p);
        }
        particles.shrink_to_fit();
    }
    for (auto& p : particles) {
        p.species = speciesNum;
    }

    return particles;
}

template<typename TParticle, typename T = typename TParticle::value_type>
std::vector<ParticleOriented<T>> initialParticlesOriented(const unsigned int& particlesNum, const std::vector<Species<T>>& allSpecies, int& speciesNum, const T& L, std::mt19937& gen, const std::string& configuration) {
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