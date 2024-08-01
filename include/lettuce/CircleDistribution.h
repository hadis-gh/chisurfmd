#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include "Circle.h"
#include "lettuce/Particle.h"
#include "lettuce/CirclesIntersectionFuncs.h"

template<typename T>
Circle<T> startCircleRandom(const T& radius, const T& areaWidth, std::mt19937 &gen) {   
    Circle<T> newCircle;
    newCircle.r = radius;

    std::uniform_real_distribution<> randomAngle(0, 2 * M_PI);
    T angle = randomAngle(gen);

    newCircle.c[0] = areaWidth * cos(angle);
    newCircle.c[1] = areaWidth * sin(angle);

    return newCircle;
}

template<typename T>
Vec<T> shootToCenter(const Circle<T> &startCircle, const T &areaWidth) {
    return -startCircle.c + areaWidth / static_cast<T>(2);
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

template<typename T>
std::vector<Particle<T>> distRandomParticles(const int& particlesNum, const T& L, const T& radius, std::mt19937& gen) {
    std::vector<Particle<T>> randomParticles(particlesNum);
    std::vector<Circle<T>> randomCircles = distRandomCircles(particlesNum, L, radius, gen);
    for (int i = 0; i < particlesNum; ++i) {
        randomParticles[i].r = randomCircles[i].c;
    }
    return randomParticles;
}

template<typename T>
std::vector<Circle<T>> distCirclesDLA(const int& shootNum, const T& L, const T& radius, std::mt19937& gen) {
    Circle<T> target = {{{L / 2, L / 2}}, radius};
    std::vector<Circle<T>> finalCircles;
    finalCircles.push_back(target);

    for (int j = 0; j < shootNum; ++j) {
        Circle<T> newCircle;
        do {
            newCircle = startCircleRandom(radius, L, gen);
        } while (has_overlap(newCircle, finalCircles));

        Vec<T> direction = shootToCenter(newCircle, L);

        Circle<T> endPoint = findStopPointAll(newCircle, direction, finalCircles);
        if (!std::isnan(endPoint.c[0])) {
            finalCircles.push_back(endPoint);
        }
    }
    return finalCircles;
}

template<typename T>
std::vector<Particle<T>> distParticleDLA(const int& shootNum, const T& L, const T& radius, std::mt19937& gen) {
    std::vector<Circle<T>> finalCircles = distCirclesDLA(shootNum, L, radius, gen);
    std::vector<Particle<T>> finalParticles(finalCircles.size());

    for (size_t i = 0; i < finalCircles.size(); ++i) {
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

template<typename T>
void writeParticle(const std::vector<Particle<T>>& particles, T radius, const std::string& fname) {
    std::ofstream output_file(fname);
    for (const auto& p : particles) {
        output_file << p.r[0] << ", " << p.r[1] << ", " << radius << "\n";
    }
}

template<typename T>
std::vector<Particle<T>> initialParticles(const unsigned int& particlesNum, const std::vector<Species<T>>& allSpecies, int& speciesNum, const T& L, std::mt19937& gen, const std::string& configuration) {
    std::vector<Particle<T>> particles;

    std::cout << "Configuration: " << configuration << std::endl;

    if (configuration == "RANDOM") {
        particles = distRandomParticles(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } else if (configuration == "DLA") {
        particles = distParticleDLA(particlesNum, L, allSpecies[speciesNum].radius, gen);
    } else if (configuration == "THREE") {
        particles.reserve(4);
        particles.push_back({{{0.0, 0.0}}, {{0.0, 0.0}}, 0});
        particles.push_back({{{1.2, 0.0}}, {{0.0, 0.0}}, 0});
        particles.push_back({{{0.0, 1.5}}, {{0.0, 0.0}}, 0});
        particles.push_back({{{1.2, 1.2}}, {{0.0, 0.0}}, 0});
    } else {
        particles.reserve(1024);
        std::ifstream config(configuration);
        std::string line;
        while (std::getline(config, line)) {
            Particle<T> p;
            std::istringstream is(line);
            is >> p.r[0] >> p.r[1];
            if (!is) {
                std::ostringstream os;
                os << "Invalid particle line " << particles.size() + 1;
                throw std::runtime_error(std::move(os).str());
            }
            // filter out particle outside of LxL box
            if (p.r[0] > L || p.r[0] < 0. || p.r[1] > L || p.r[1] < 0.) {
                continue;
            }
            is >> p.v[0] >> p.v[1];
            if (!is) {
                p.v = {{0., 0.}};
            }
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