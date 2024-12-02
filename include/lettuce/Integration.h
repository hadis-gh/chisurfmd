#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/LennardJones.h"
#include "lettuce/LennardJonesOriented.h"
#include "lettuce/Thermostat.h"

template<typename TParticle, typename Force, typename T = typename TParticle::value_type>
void EulerSymplecticStep(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto r = getGeneralizedPositions(particles[i]);
        const auto v = getGeneralizedVelocities(particles[i]);
        setGeneralizedVelocities(particles[i], v + accelerations[i] * dt);
        setGeneralizedPositions(particles[i], r + particles[i].v * dt);
        implementPBC(particles[i], boxPBC);
    }
}

template<typename TParticle, typename Force, typename T = typename TParticle::value_type>
void EulerStep(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto r = getGeneralizedPositions(particles[i]);
        const auto v = getGeneralizedVelocities(particles[i]);
        setGeneralizedPositions(particles[i], r + v * dt);
        setGeneralizedVelocities(particles[i], v + accelerations[i] * dt);
        implementPBC(particles[i], boxPBC);
    }
}

template<typename TParticle, typename Force, typename T = typename TParticle::value_type>
void VelocityVerletStep(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto old_accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto r = getGeneralizedPositions(particles[i]);
        const auto v = getGeneralizedVelocities(particles[i]);

        setGeneralizedPositions(particles[i], r + v * dt + old_accelerations[i] * dt * dt / 2.0);
        implementPBC(particles[i], boxPBC);
    }

    const auto new_accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto v = getGeneralizedVelocities(particles[i]);
        setGeneralizedVelocities(particles[i], v + (old_accelerations[i] + new_accelerations[i]) * dt / 2.0);
    }
}

template<typename TParticle, typename Integrator, typename Force, typename T = typename TParticle::value_type>
void integrate(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, T dt, T Time, const T& boxPBC, Force&& force, Integrator&& integrator) {
    int numSteps = static_cast<int>(Time / dt);

    for (int i = 0; i < numSteps; ++i) {
        integrator(particles, allSpecies, dt, boxPBC, std::forward<Force>(force));
    }
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calAverageNeighbors(const std::vector<TParticle>& particles, const T& distance) {
    T totalNeighbors = 0;
    for (const auto& p1 : particles) {
        int eachParticleNeighbors = 0;
        for (const auto& p2 : particles) {
            if (p1.r != p2.r && (p1.r - p2.r).abs() <= distance) {
                ++eachParticleNeighbors;
            }
        }
        totalNeighbors += eachParticleNeighbors;
    }
    return totalNeighbors / particles.size();
}

template<typename TParticle, typename T = typename TParticle::value_type>
Vec<T> calCOMposition(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comPos;
    T totalMass = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        comPos += p.r * mass;
        totalMass += mass;
    }
    return comPos / totalMass;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calAngularMomentum2D(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const Vec<T>& comPos) {
    T L = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        Vec<T> r_com = p.r - comPos;
        L += mass * (r_com[0] * p.v[1] - r_com[1] * p.v[0]);
    }
    return L;
}

template<typename TParticle, typename T = typename TParticle::value_type>
T calMomentOfInertia2D(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const Vec<T>& comPos) {
    T I = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        Vec<T> r_com = p.r - comPos;
        I += mass * r_com.abs2();
    }
    return I;
}

template<typename TParticle, typename T = typename TParticle::value_type>
void removeCOMvelocityRotation2D(std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comPos = calCOMposition(particles, allSpecies);
    T angularMomentum = calAngularMomentum2D(particles, allSpecies, comPos);
    T momentOfInertia = calMomentOfInertia2D(particles, allSpecies, comPos);

    T angularVelocity = angularMomentum / momentOfInertia;
    // there is possiblity to remove mass from these eqs - try to simplify it

    for (auto& p : particles) {
        Vec<T> r_com = p.r - comPos;
        Vec<T> v_rot = {(-r_com[1] * angularVelocity, r_com[0] * angularVelocity)};
        p.v -= v_rot;
    }
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeInitialParticles(const std::vector<TParticle>& particles, const T radius) {
    std::ofstream initialParticles("FirstConfigPlot.dat");
    for (const auto& p : particles) {
        initialParticles << p.r[0] << " " << p.r[1] << " " << radius << std::endl;
    }
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writePositionToFile(const std::vector<TParticle>& particles, std::ostream& file, const T& dt, const int& step) {
    file << dt * step << " ";
    for (const auto& p : particles) {
        // constexpr int D = degreesOfFreedom<TParticle>();
        auto q = getGeneralizedPositions(p);
        for (int a = 0; a < q.size(); ++a) {
            file << p.r[a] << " ";
        }
    }
    file << "\n";
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeKineticEToFile(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file) {
    for (const auto& p : particles) {
        // auto q = getGeneralizedPositions(p);
        // for (int a = 0; a < q.size() - 1; ++a) {
        file << calParticleKineticEnergy(p, allSpecies) << " ";
        // }
    }
    file << "\n";
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeRelativeKineticEToFile(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file) {
    Vec<T> comVelocity = comVelocity(particles, allSpecies);

    for (const auto& p : particles) {
        Vec<T> relativeVelocity = p.v - comVelocity;
        file << calParticleRelativeKineticEnergy(p, allSpecies, relativeVelocity) << " ";
    }
    file << "\n";
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file) {
    file << calInternalTemperature(particles, allSpecies) << std::endl;
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeRealTemperature(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file) {
    file << calRawTemperature(particles, allSpecies) << std::endl;
}

template<typename TParticle, typename Potential, typename T = typename TParticle::value_type>
void writePotentialEToFile(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Potential&& potential, std::ostream& file) {
    const auto potentialEnergy = calAllAccelerations(particles, allSpecies, boxPBC, std::forward<Potential>(potential));
    for (const auto& u : potentialEnergy) {
        file << u.abs() << " ";
    }
    file << "\n";
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeAverageNeighborToFile(const std::vector<TParticle>& particles, const std::vector<T>& distances, std::ostream& file, const T& dt, const int& step) {
    file << dt * step << " ";
    for (const T& d : distances) {
        file << calAverageNeighbors(particles, d) << " ";
    }
    file << "\n";
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeAverageNeighborToFile(const std::vector<TParticle>& particles, const T& distance, std::ostream& file, const T& dt, const int& step) {
    file << dt * step << " " << calAverageNeighbors(particles, distance) << "\n";
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeComVelocityToFile(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file, const T& dt, const int& step) {
    file << dt * step << " " << calCOMVelocity(particles, allSpecies) << "\n";
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeComAngularVelocityToFile(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file, const T& dt, const int& step) {
    Vec<T> comPos = calCOMposition(particles, allSpecies);

    file << dt * step << " " << calAngularMomentum2D(particles, allSpecies, comPos) << "\n";
}