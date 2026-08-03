#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <memory>
#include <tuple>

#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/ParticleDot.h"
#include "chisurfmd/core/ParticleOriented.h"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

template<typename T>
using ReferencePotentialList =
    std::vector<std::pair<T, T>>;  // (r, energy)


enum class RefInteractionType {
    EP = 0,
    OP = 1,
    EA = 2,
    OA = 3
};


// ==================================
// Helper Functions: TabularDFT
// ==================================

template<typename T>
void readReferencePotentialFile(
    const std::string& potentialPath,
    ReferencePotentialList<T>& potential,
    RefInteractionType interactionType,
    T distanceScale,
    T energyScale
) {
    potential.clear();

    T phi1;
    T phi2;
    T dz;
    T r;
    T energy;

    int count = 0;

    const auto interactionIndex =
        static_cast<int>(interactionType);

    std::ifstream fileStream(
        potentialPath
        + std::to_string(interactionIndex)
        + ".dat"
    );

    while (fileStream >> phi1 >> phi2 >> dz >> r >> energy) {
        if (
            std::abs(dz) > 1e-6
            || std::abs(phi1) > 1e-6
            || std::abs(phi2) > 1e-6
        ) {
            continue;
        }

        const T scaledDistance =
            r * distanceScale;

        potential.push_back({
            scaledDistance,
            energy
        });

        count++;
    }

    if (count == 0) {
        for (T distance = 7.9; distance <= 20.0; distance += 0.1) {
            potential.push_back({
                distance * distanceScale,
                0.0
            });
        }
    }

    std::sort(
        potential.begin(),
        potential.end(),
        [](const auto& first, const auto& second) {
            return first.first < second.first;
        }
    );
}


template<typename T>
T interpolateReferenceEnergy(
    const ReferencePotentialList<T>& potential,
    const T& r,
    const T& energyScale
) {
    if (potential.empty()) {
        return T{0};
    }

    if (r < potential.front().first) {
        const T minimumDistance =
            potential.front().first;

        return
            T{1000}
            * energyScale
            * std::exp(
                T{10} * (minimumDistance - r)
            );
    }

    if (r > potential.back().first) {
        return T{0};
    }

    for (
        std::size_t index = 0;
        index + 1 < potential.size();
        ++index
    ) {
        if (
            r >= potential[index].first
            && r <= potential[index + 1].first
        ) {
            const T r1 = potential[index].first;
            const T e1 = potential[index].second;

            const T r2 = potential[index + 1].first;
            const T e2 = potential[index + 1].second;

            if (std::abs(r2 - r1) < T{1e-10}) {
                return e1 * energyScale;
            }

            const T interpolationFraction =
                (r - r1) / (r2 - r1);

            const T energy =
                e1 * (T{1} - interpolationFraction)
                + e2 * interpolationFraction;

            return energy * energyScale;
        }
    }

    return potential.back().second * energyScale;
}


template<typename T>
T interpolateReferenceDerivative(
    const ReferencePotentialList<T>& potential,
    const T& r,
    const T& energyScale
) {
    if (potential.empty()) {
        return T{0};
    }

    const T epsilon = T{1e-6};

    const T lowerEnergy =
        interpolateReferenceEnergy(
            potential,
            r - epsilon,
            energyScale
        );

    const T upperEnergy =
        interpolateReferenceEnergy(
            potential,
            r + epsilon,
            energyScale
        );

    return
        (upperEnergy - lowerEnergy)
        / (T{2} * epsilon);
}


template<typename T>
RefInteractionType getRefInteractionType(
    const ParticleOriented<T>& p1,
    const ParticleOriented<T>& p2
) {
    const bool equalHandedness =
        p1.handedness * p2.handedness > 0;

    const bool parallelAlignment =
        p1.alignment * p2.alignment > 0;

    if (equalHandedness && parallelAlignment) {
        return RefInteractionType::EP;
    }

    if (!equalHandedness && parallelAlignment) {
        return RefInteractionType::OP;
    }

    if (equalHandedness && !parallelAlignment) {
        return RefInteractionType::EA;
    }

    return RefInteractionType::OA;
}


// ==================================
// Potential Calculation Class: TabularDFT
// ==================================

template<typename T>
class TabularDFTPotential {
public:
    TabularDFTPotential(
        const std::string& potentialPath,
        T energyScale,
        T distanceScale,
        T cutoff
    )
        : m_potentialPath(potentialPath),
          m_energyScale(energyScale),
          m_distanceScale(distanceScale),
          m_cutoff(cutoff)
    {
        readReferencePotentialFile(
            m_potentialPath,
            m_referenceEP,
            RefInteractionType::EP,
            m_distanceScale,
            m_energyScale
        );

        readReferencePotentialFile(
            m_potentialPath,
            m_referenceOP,
            RefInteractionType::OP,
            m_distanceScale,
            m_energyScale
        );

        readReferencePotentialFile(
            m_potentialPath,
            m_referenceEA,
            RefInteractionType::EA,
            m_distanceScale,
            m_energyScale
        );

        readReferencePotentialFile(
            m_potentialPath,
            m_referenceOA,
            RefInteractionType::OA,
            m_distanceScale,
            m_energyScale
        );
    }


    T operator()(
        const ParticleOriented<T>& p1,
        const ParticleOriented<T>& p2,
        const Vec<T, 2>&,
        const T R
    ) const {
        const T scaledDistance =
            R * m_distanceScale;

        if (scaledDistance > m_cutoff || R <= T{0}) {
            return T{0};
        }

        const auto& potential =
            getReferencePotential(p1, p2);

        if (potential.empty()) {
            return T{0};
        }

        return interpolateReferenceEnergy(
            potential,
            scaledDistance,
            m_energyScale
        );
    }


    std::pair<T, Vec<T, 3>> energyAndForces(
        const ParticleOriented<T>& p1,
        const ParticleOriented<T>& p2,
        const Vec<T, 2>& dr,
        const T R
    ) const {
        const T scaledDistance =
            R * m_distanceScale;

        if (scaledDistance > m_cutoff || R <= T{0}) {
            return {
                T{0},
                {{T{0}, T{0}, T{0}}}
            };
        }

        const auto& potential =
            getReferencePotential(p1, p2);

        if (potential.empty()) {
            return {
                T{0},
                {{T{0}, T{0}, T{0}}}
            };
        }

        const T energy =
            interpolateReferenceEnergy(
                potential,
                scaledDistance,
                m_energyScale
            );

        const T energyDerivative =
            interpolateReferenceDerivative(
                potential,
                scaledDistance,
                m_energyScale
            );

        T forceX = T{0};
        T forceY = T{0};

        if (R > T{1e-10}) {
            forceX =
                -energyDerivative * dr[0] / R;

            forceY =
                -energyDerivative * dr[1] / R;
        }

        const T torque = T{0};

        return {
            energy,
            {{forceX, forceY, torque}}
        };
    }


private:
    const ReferencePotentialList<T>& getReferencePotential(
        const ParticleOriented<T>& p1,
        const ParticleOriented<T>& p2
    ) const {
        const RefInteractionType type =
            getRefInteractionType(p1, p2);

        switch (type) {
            case RefInteractionType::EP:
                return m_referenceEP;

            case RefInteractionType::OP:
                return m_referenceOP;

            case RefInteractionType::EA:
                return m_referenceEA;

            case RefInteractionType::OA:
                return m_referenceOA;
        }

        throw std::runtime_error(
            "Unknown reference interaction type"
        );
    }


    std::string m_potentialPath;

    T m_energyScale;
    T m_distanceScale;
    T m_cutoff;

    ReferencePotentialList<T> m_referenceEP;
    ReferencePotentialList<T> m_referenceOP;
    ReferencePotentialList<T> m_referenceEA;
    ReferencePotentialList<T> m_referenceOA;
};


// ==================================
// Force Calculation Class: TabularDFT
// ==================================

template<
    typename TParticle,
    typename T = typename TParticle::value_type
>
class TabularDFTForce {
public:
    using value_type = T;

    TabularDFTForce(
        const std::string& potentialPath,
        T energyScale,
        T distanceScale,
        T cutoff
    )
        : m_potential(
              potentialPath,
              energyScale,
              distanceScale,
              cutoff
          )
    {
    }


    Vec<T, 3> operator()(
        const TParticle& p1,
        const TParticle& p2,
        const Vec<T, 2>& dr,
        const T r
    ) const {
        return m_potential
            .energyAndForces(p1, p2, dr, r)
            .second;
    }


private:
    TabularDFTPotential<T> m_potential;
};


// ==================================
// Factory Struct for TabularDFT
// ==================================

template<typename Particle, typename SFINAE = void>
struct TabularDFT;


template<typename T>
struct TabularDFT<ParticleOriented<T>> {
    static void initProgramOptions(
        po::options_description& desc
    ) {
        IsotropicLJ<ParticleDot<T>>::initProgramOptions(desc);

        desc.add_options()
            (
                "dftb.potentialPath",
                po::value<std::string>()->default_value(
                    "/home/hadis/vector/DTFB_files/"
                ),
                "DFTB potential path"
            )
            (
                "dftb.energyScale",
                po::value<T>()->default_value(0.0001),
                "Energy scaling factor"
            )
            (
                "dftb.distanceScale",
                po::value<T>()->default_value(1.0),
                "Distance scaling factor"
            );
    }


    static auto force(
        const po::variables_map& vm
    ) {
        try {
            return TabularDFTForce<ParticleOriented<T>>(
                vm["dftb.potentialPath"].as<std::string>(),
                vm["dftb.energyScale"].as<T>(),
                vm["dftb.distanceScale"].as<T>(),
                vm["LJcutoff"].as<T>()
            );
        } catch (const boost::bad_any_cast& error) {
            std::cerr
                << "Error initializing TabularDFTForce: "
                << error.what()
                << std::endl;

            throw;
        }
    }


    static auto potential(
        const po::variables_map& vm
    ) {
        try {
            return TabularDFTPotential<T>(
                vm["dftb.potentialPath"].as<std::string>(),
                vm["dftb.energyScale"].as<T>(),
                vm["dftb.distanceScale"].as<T>(),
                vm["LJcutoff"].as<T>()
            );
        } catch (const boost::bad_any_cast& error) {
            std::cerr
                << "Error initializing TabularDFTPotential: "
                << error.what()
                << std::endl;

            throw;
        }
    }


    using ForceType =
        TabularDFTForce<ParticleOriented<T>>;
};