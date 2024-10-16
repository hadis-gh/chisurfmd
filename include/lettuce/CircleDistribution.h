#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <sstream>
#include <stdexcept>
#include "Circle.h"
#include "lettuce/Particle.h"
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
std::vector<TParticle> manualRandomParticles(const int& particlesNum, const T& L, const T& radius, std::mt19937& gen) {
    int topSquareRoot = static_cast<int>(std::ceil(std::sqrt(particlesNum)));
    std::vector<TParticle> particles;
    T currDensity = packingDensity(topSquareRoot * topSquareRoot, radius, L);
    const double maxPackingDensity = 0.7854;

    int maxParticles = static_cast<int>(maxPackingDensity * L * L / (M_PI * radius * radius));

    if (currDensity > maxPackingDensity) {
        std::cout << "maximum possible number for particles is : " << maxParticles << " but you entered: " << particlesNum << std::endl;
        throw std::runtime_error("Too many particles for this area!");
    }

    particles.reserve(particlesNum);

    T distance = L / (topSquareRoot + 1);

    const T dr = distance / 2 - radius;

    for (int i = 0; i < topSquareRoot; ++i) {
        for (int j = 0; j < topSquareRoot; ++j) {
            if (particles.size() < particlesNum) {
                auto newParticle = createRandomParticle<TParticle>(gen);
                for(int a = 0; a < 2; ++a)
                {
                    newParticle.r[a] = ((newParticle.r[a] * 2) - 1) * dr;
                    newParticle.r[a] += (i + 1) * distance;
                }
                particles.push_back(newParticle);
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
        // particles = distRandomParticles(particlesNum, L, allSpecies[speciesNum].radius, gen);
        particles = manualRandomParticles<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } else if (configuration == "DLA") {
        particles = distParticleDLA<TParticle>(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } 
    else {
        particles.reserve(1024);
        std::ifstream config(configuration);
        std::string line;
        while (std::getline(config, line)) {
            TParticle p;
            std::istringstream is(line);
            constexpr int D = degreesOfFreedom<TParticle>();
            Vec<T, D> q;
            for(int a = 0; a < D; ++a)
            {
                is >> q[a];
            }
            setGeneralizedPositions(p, q);
            if (!is) {
                std::ostringstream os;
                os << "Invalid particle line " << particles.size() + 1;
                throw std::runtime_error(std::move(os).str());
            }
            // filter out particle outside of LxL box
            if (p.r[0] > L || p.r[0] < 0. || p.r[1] > L || p.r[1] < 0.) {
                continue;
            }
            for(int a = 0; a < D; ++a)
            {
                is >> q[a];
            }
            if (!is) {
                q.fill(0);
            }
            setGeneralizedVelocities(p, q);
            if (!is.eof()) {
                std::ostringstream os;
                os << "Invalid particle line " << particles.size() + 1 << ": unread characters.";
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