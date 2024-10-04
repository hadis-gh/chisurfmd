#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/LennardJones.h"
#include "lettuce/Thermostat.h"

template<typename T>
void implementPBC(Particle<T>& p, const T& boxPBC) {
    for (int i = 0; i < 2; ++i) {
        if (p.r[i] > boxPBC) { p.r[i] -= boxPBC; }
        else if (p.r[i] < 0) { p.r[i] += boxPBC; }
    }
}

template<typename T, typename Force>
void EulerSymplecticStep(std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].v += accelerations[i] * dt;
        particles[i].r += particles[i].v * dt;
        implementPBC(particles[i], boxPBC);
    }
}

template<typename T, typename Force>
void EulerStep(std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].r += particles[i].v * dt;
        particles[i].v += accelerations[i] * dt;
        implementPBC(particles[i], boxPBC);
    }
}

template<typename T, typename Force>
void VelocityVerletStep(std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, const T& dt, const T& boxPBC, const Force& force) {
    const auto old_accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);

    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].r += particles[i].v * dt + old_accelerations[i] * dt * dt / 2.0;
        implementPBC(particles[i], boxPBC);
    }

    const auto new_accelerations = calAllAccelerations(particles, allSpecies, boxPBC, force);
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].v += (old_accelerations[i] + new_accelerations[i]) * dt / 2.0;
    }
}

template<typename T, typename Integrator, typename Force>
void integrate(std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, T dt, T Time, const T& boxPBC, Force&& force, Integrator&& integrator) {
    int numSteps = static_cast<int>(Time / dt);

    for (int i = 0; i < numSteps; ++i) {
        integrator(particles, allSpecies, dt, boxPBC, std::forward<Force>(force));
    }
}

template<typename T>
T calAverageNeighbors(const std::vector<Particle<T>>& particles, const T& distance) {
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

template<typename T>
Vec<T> calCOMposition(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comPos;
    T totalMass = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        comPos += p.r * mass;
        totalMass += mass;
    }
    return comPos / totalMass;
}

template<typename T>
Vec<T> calAngularMomentum(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, const Vec<T>& comPos) {
    Vec<T> L;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        Vec<T> r_com = p.r - comPos;
        L += r_com.cross(p.v) * mass;
    }
    return L;
}

template<typename T>
T calMomentOfInertia(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, const Vec<T>& com) {
    T I = 0;
    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        Vec<T> r_com = p.r - com;
        I += mass * r_com.abs2();
    }
    return I;
}

template<typename T>
void removeCOMvelocityRotation(std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies) {
    Vec<T> comPos = calCOMposition(particles, allSpecies);
    Vec<T> angularMomentum = calAngularMomentum(particles, allSpecies, comPos);
    T momentOfInertia = calMomentOfInertia(particles, allSpecies, comPos);

    Vec<T> angularVelocity = angularMomentum / momentOfInertia;  // Angular velocity ω = L / I

    for (auto& p : particles) {
        Vec<T> r_com = p.r - comPos;
        Vec<T> v_rot = angularVelocity.cross(r_com);  // Rotational velocity component
        p.v -= v_rot;
    }
}

template<typename T>
void writeInitialParticles(const std::vector<Particle<T>>& particles, const T radius) {
    std::ofstream initialParticles("FirstConfigPlot.dat");
    for (const auto& p : particles) {
        initialParticles << p.r[0] << " " << p.r[1] << " " << radius << std::endl;
    }
}

template<typename T>
void writePositionToFile(const std::vector<Particle<T>>& particles, std::ostream& file, const T& dt, const int& step) {
    file << dt * step << " ";
    for (const auto& p : particles) {
        file << p.r[0] << " " << p.r[1] << " ";
    }
    file << "\n";
}

template<typename T>
void writeKineticEToFile(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file) {
    for (const auto& p : particles) {
        file << calParticleKineticEnergy(p, allSpecies) << " ";
    }
    file << "\n";
}

template<typename T>
void writeRelativeKineticEToFile(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file) {
    Vec<T> comVelocity = comVelocity(particles, allSpecies);

    for (const auto& p : particles) {
        Vec<T> relativeVelocity = p.v - comVelocity;
        file << calParticleRelativeKineticEnergy(p, allSpecies, relativeVelocity) << " ";
    }
    file << "\n";
}

template<typename T>
void writeTemperature(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file) {
    file << calInternalTemperature(particles, allSpecies) << std::endl;
}

template<typename T>
void writeRealTemperature(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file) {
    file << calRawTemperature(particles, allSpecies) << std::endl;
}

template<typename T, typename Potential>
void writePotentialEToFile(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, const T& boxPBC, Potential&& potential, std::ostream& file) {
    const auto potentialEnergy = calAllAccelerations(particles, allSpecies, boxPBC, std::forward<Potential>(potential));
    for (const auto& u : potentialEnergy) {
        file << u.abs() << " ";
    }
    file << "\n";
}

template<typename T>
void writeAverageNeighborToFile(const std::vector<Particle<T>>& particles, const std::vector<T>& distances, std::ostream& file, const T& dt, const int& step) {
    file << dt * step << " ";
    for (const T& d : distances) {
        file << calAverageNeighbors(particles, d) << " ";
    }
    file << "\n";
}

template<typename T>
void writeAverageNeighborToFile(const std::vector<Particle<T>>& particles, const T& distance, std::ostream& file, const T& dt, const int& step) {
    file << dt * step << " " << calAverageNeighbors(particles, distance) << "\n";
}

template<typename T>
void writeComVelocityToFile(const std::vector<Particle<T>>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file, const T& dt, const int& step) {
    file << dt * step << " " << calCOMVelocity(particles, allSpecies) << "\n";
}
