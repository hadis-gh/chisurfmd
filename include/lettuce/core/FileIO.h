#pragma once

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "lettuce/core/Vec.h"

template<typename TParticle, typename T = typename TParticle::value_type>
void writeInitialParticles(const std::vector<TParticle>& particles, const T radius) {
    std::ofstream initialParticles("FirstConfigPlot.dat");
    for (const auto& p : particles) {
        auto q = getGeneralizedPositions(p);
        for (int a = 0; a < q.size(); ++ a){
            initialParticles << q[a] << " ";
        }
        initialParticles << std::endl;
    }
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writePositionToFile(const std::vector<TParticle>& particles, std::ostream& file, const T& dt, const int& step) {
    file << dt * step << " ";
    for (const auto& p : particles) {
        auto q = getGeneralizedPositions(p);
        for (int a = 0; a < q.size(); ++a) {
            file << q[a] << " ";
        }
    }
    file << "\n";
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeKineticEToFile(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, adios2::fstream& oStream) {
	std::vector<T> eKin(particles.size());
    for (int a = 0; a < particles.size(); ++a) {
        eKin[a] = calParticleKineticEnergy(particles[a], allSpecies);
    }
	oStream.write("kineticEnergy", eKin.data(), {particles.size()},
			{std::size_t(0u)}, {particles.size()});
}

template<typename TParticle, typename T = typename TParticle::value_type>
void writeAverageKineticEnergies(const std::vector<TParticle>& particles, const std::vector<Species<T>>& allSpecies, std::ostream& file) {
    T totalKE_x = 0.0;
    T totalKE_y = 0.0;
    T totalKE_phi = 0.0;

    size_t numParticles = particles.size();

    for (const auto& p : particles) {
        const T mass = allSpecies[p.species].mass;
        const T MOI = allSpecies[p.species].momentOfInertia;

        auto vel = getGeneralizedVelocities(p);

        totalKE_x += 0.5 * mass * vel[0] * vel[0];
        totalKE_y += 0.5 * mass * vel[1] * vel[1];

        if (vel.size() > 2) {
            totalKE_phi += 0.5 * MOI * vel[2] * vel[2];
        }
    }

    T avgKE_x = totalKE_x / numParticles;
    T avgKE_y = totalKE_y / numParticles;
    T avgKE_phi = (totalKE_phi > 0) ? (totalKE_phi / numParticles) : 0.0;

    file << avgKE_x << " " << avgKE_y << " " << avgKE_phi << "\n";
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
    for (auto &p1 : particles) {
        double pot = 0.;
        for (auto &p : particles) {
            if (p1.r != p.r) {
                Vec<T> dr = p.r - p1.r;

                for (int i = 0; i < 2; ++i) {
                    if (dr[i] > boxPBC / 2) { dr[i] -= boxPBC; }
                    else if (dr[i] < -boxPBC / 2) { dr[i] += boxPBC; }
                }

                const T r = dr.abs();
                if (r == 0) continue;

                if constexpr (std::is_same_v<TParticle, ParticleOriented<T>>) {
                    T deltaPhi = p.phi - p1.phi;
                    pot += potential(r, deltaPhi);
                } else {
                    pot += potential(r);
                }
            }
        }
        file << pot << "\t";
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

template<typename TParticle, typename T = typename TParticle::value_type>
void saveParticlesWithVelocities(std::ostream& output, const std::vector<TParticle>& particles) {
    for (const auto& p : particles) {
        auto q = getGeneralizedPositions(p);
        auto qDot = getGeneralizedVelocities(p);
        constexpr int D = degreesOfFreedom<TParticle>();

        for (int a = 0; a < D; ++a) {
            output << std::setprecision(13) << std::scientific << q[a] << "\t";
        }
        for (int a = 0; a < D - 1; ++a) {
            output << std::setprecision(13) << std::scientific << qDot[a] << "\t";
        }
        output << std::setprecision(13) << std::scientific << qDot[D - 1] << std::endl;
    }
}