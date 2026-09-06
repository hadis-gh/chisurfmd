#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <boost/program_options.hpp>
#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/ParticleOriented.h"
namespace po = boost::program_options;
// ============================================================
// Interaction types and configuration
// ============================================================
enum class RefInteractionType { EP, EA, OP, OA };
inline const char* interactionName(RefInteractionType type) {
    switch (type) {
        case RefInteractionType::EP: return "EP";
        case RefInteractionType::EA: return "EA";
        case RefInteractionType::OP: return "OP";
        case RefInteractionType::OA: return "OA";
    }
    throw std::runtime_error("Unknown TabularDFT interaction type");
}
inline std::string interactionFilename(RefInteractionType type) {
    return std::string("E_all_") + interactionName(type) + ".dat";
}
template<typename T>
RefInteractionType getRefInteractionType(
    const ParticleOriented<T>& p1,
    const ParticleOriented<T>& p2
) {
    const bool equal = p1.handedness == p2.handedness;
    const bool parallel = p1.alignment == p2.alignment;
    if (equal) return parallel ? RefInteractionType::EP : RefInteractionType::EA;
    return parallel ? RefInteractionType::OP : RefInteractionType::OA;
}
template<typename T>
struct TabularDFTConfig {
    T isolatedHelixEnergyEV;
    bool inputEnergyIsBinding;
    T E0;                    // eV
    T L0;                    // Angstrom
    T cutoff;                // reduced distance
    T zeta;                  // Angstrom
    T zetaTolerance;
    int screwStepDeg;
    bool requireFullAngularCoverage;
};

// ------------------------------------------------------------
// Fixed reference-data assumptions.
// Change them here if a future dataset uses different units/settings.
// ------------------------------------------------------------
namespace TabularDFTDefaults {
inline constexpr long double isolatedHelixEnergyEV = -6227.1749L;
inline constexpr bool inputEnergyIsBinding = false;
inline constexpr long double E0 = 0.4L;               // eV
inline constexpr long double L0 = 7.0L;               // Angstrom
inline constexpr long double cutoffAngstrom = 20.0L;  // Angstrom
inline constexpr long double zeta = 0.0L;             // Angstrom
inline constexpr long double zetaTolerance = 1e-8L;
inline constexpr int screwStepDeg = 20;
inline constexpr bool requireFullAngularCoverage = true;
}

// ------------------------------------------------------------
inline int inferScrewDirection(RefInteractionType type) {
    // Equal-handed interactions:    EP, EA -> +1
    // Opposite-handed interactions: OP, OA -> -1
    switch (type) {
        case RefInteractionType::EP:
        case RefInteractionType::EA:
            return +1;

        case RefInteractionType::OP:
        case RefInteractionType::OA:
            return -1;
    }

    throw std::runtime_error("Unknown TabularDFT interaction type");
}

// ------------------------------------------------------------
template<typename T>
TabularDFTConfig<T> defaultTabularDFTConfig() {
    return {
        static_cast<T>(TabularDFTDefaults::isolatedHelixEnergyEV),
        TabularDFTDefaults::inputEnergyIsBinding,
        static_cast<T>(TabularDFTDefaults::E0),
        static_cast<T>(TabularDFTDefaults::L0),
        static_cast<T>(
            TabularDFTDefaults::cutoffAngstrom /
            TabularDFTDefaults::L0
        ),
        static_cast<T>(TabularDFTDefaults::zeta),
        static_cast<T>(TabularDFTDefaults::zetaTolerance),
        TabularDFTDefaults::screwStepDeg,
        TabularDFTDefaults::requireFullAngularCoverage
    };
}
// ============================================================
// Small helpers and data types
// ============================================================
inline int positiveModulo(int x, int period) {
    const int r = x % period;
    return r < 0 ? r + period : r;
}
template<typename T>
T radiansToWrappedDegrees(T angle) {
    constexpr T degPerRad =
        static_cast<T>(180.0L / 3.141592653589793238462643383279502884L);
    angle = std::fmod(angle * degPerRad, T{360});
    return angle < T{0} ? angle + T{360} : angle;
}
template<typename T>
struct RadialPoint {
    T r;       // R / L0
    T energy;  // E_binding / E0
};
template<typename T>
using RadialCurve = std::vector<RadialPoint<T>>;
struct AngularKey {
    int phi1 = 0;
    int phi2 = 0;
    bool operator<(const AngularKey& other) const {
        return phi1 != other.phi1 ? phi1 < other.phi1 : phi2 < other.phi2;
    }
    bool operator==(const AngularKey& other) const {
        return phi1 == other.phi1 && phi2 == other.phi2;
    }
};
template<typename T>
struct TabulatedValue {
    T energy = T{0};
    T dUdr = T{0};
    T dUdPhi1 = T{0}; // per radian
    T dUdPhi2 = T{0}; // per radian
};
// ============================================================
// One screw-periodic reference table
// ============================================================
template<typename T>
class ScrewPeriodicReferenceTable {
public:
    ScrewPeriodicReferenceTable(
        const std::string& filename,
        const TabularDFTConfig<T>& cfg,
        int screwDirection
    )
        : m_expandedToRaw(360 * 360),
          m_screwStepDeg(cfg.screwStepDeg),
          m_screwDirection(screwDirection)
    {
        validate(cfg);
        readRawTable(filename, cfg);
        buildScrewCoverage();
        if (cfg.requireFullAngularCoverage && m_missingAngularNodes > 0) {
            throw std::runtime_error(
                "Incomplete screw-expanded angular grid in " + filename +
                ": missing " + std::to_string(m_missingAngularNodes) +
                " of 129600 nodes."
            );
        }
    }
    std::size_t missingAngularNodes() const {
        return m_missingAngularNodes;
    }
    // --------------------------------------------------------
    // U(phi1, phi2, r):
    //   1. convert MD angles from radians to [0,360) degrees
    //   2. evaluate the four neighboring integer angular nodes
    //   3. bilinearly interpolate energy and radial derivative
    //   4. differentiate the same bilinear surface for torques
    // --------------------------------------------------------
    TabulatedValue<T> evaluate(T r, T phi1Rad, T phi2Rad) const {
        const T phi1 = radiansToWrappedDegrees(phi1Rad);
        const T phi2 = radiansToWrappedDegrees(phi2Rad);
        const T f1 = std::floor(phi1);
        const T f2 = std::floor(phi2);
        const int p10 = static_cast<int>(f1) % 360;
        const int p20 = static_cast<int>(f2) % 360;
        const int p11 = (p10 + 1) % 360;
        const int p21 = (p20 + 1) % 360;
        const T a = phi1 - f1;
        const T b = phi2 - f2;
        const auto v00 = evaluateAngularNode(p10, p20, r);
        const auto v10 = evaluateAngularNode(p11, p20, r);
        const auto v01 = evaluateAngularNode(p10, p21, r);
        const auto v11 = evaluateAngularNode(p11, p21, r);
        TabulatedValue<T> out;
        out.energy = bilinear(
            v00.first, v10.first, v01.first, v11.first, a, b
        );
        out.dUdr = bilinear(
            v00.second, v10.second, v01.second, v11.second, a, b
        );
        const T dPhi1Deg =
            (T{1} - b) * (v10.first - v00.first) +
            b * (v11.first - v01.first);
        const T dPhi2Deg =
            (T{1} - a) * (v01.first - v00.first) +
            a * (v11.first - v10.first);
        constexpr T degPerRad =
            static_cast<T>(180.0L / 3.141592653589793238462643383279502884L);
        out.dUdPhi1 = dPhi1Deg * degPerRad;
        out.dUdPhi2 = dPhi2Deg * degPerRad;
        return out;
    }
private:
    using CurveMap = std::map<AngularKey, RadialCurve<T>>;
    CurveMap m_rawCurves;
    // Heap-backed 360 x 360 lookup. Each expanded node stores the raw
    // angular nodes that contribute to it after screw expansion.
    std::vector<std::vector<AngularKey>> m_expandedToRaw;
    int m_screwStepDeg;
    int m_screwDirection;
    std::size_t m_missingAngularNodes = 0;
    static std::size_t angularIndex(int phi1, int phi2) {
        return static_cast<std::size_t>(phi1 * 360 + phi2);
    }
    static T bilinear(
        T v00, T v10, T v01, T v11,
        T a, T b
    ) {
        return
            (T{1} - a) * (T{1} - b) * v00 +
            a * (T{1} - b) * v10 +
            (T{1} - a) * b * v01 +
            a * b * v11;
    }
    void validate(const TabularDFTConfig<T>& cfg) const {
        if (!(cfg.E0 > T{0})) {
            throw std::runtime_error("dftb.epsilonEV must be > 0");
        }
        if (!(cfg.L0 > T{0})) {
            throw std::runtime_error("dftb.sigmaAngstrom must be > 0");
        }
        if (cfg.screwStepDeg <= 0 || 360 % cfg.screwStepDeg != 0) {
            throw std::runtime_error(
                "dftb.screwStepDeg must be positive and divide 360"
            );
        }
        if (m_screwDirection != 1 && m_screwDirection != -1) {
            throw std::runtime_error("screwDirection must be +1 or -1");
        }
    }
    // --------------------------------------------------------
    // Read only the selected zeta slice and convert to reduced units:
    //
    // E_binding = E_pair - 2 E_isolated
    // r*        = R / L0
    // U*        = E_binding / E0
    // --------------------------------------------------------
    void readRawTable(
        const std::string& filename,
        const TabularDFTConfig<T>& cfg
    ) {
        std::ifstream input(filename);
        if (!input) {
            throw std::runtime_error(
                "Could not open reference potential file: " + filename
            );
        }
        T phi1, phi2, zeta, R, inputEnergy;
        std::size_t accepted = 0;
        while (input >> phi1 >> phi2 >> zeta >> R >> inputEnergy) {
            if (std::abs(zeta - cfg.zeta) > cfg.zetaTolerance) continue;
            const int p1 = static_cast<int>(std::lround(phi1));
            const int p2 = static_cast<int>(std::lround(phi2));
            if (
                std::abs(phi1 - static_cast<T>(p1)) > T{1e-6} ||
                std::abs(phi2 - static_cast<T>(p2)) > T{1e-6}
            ) {
                throw std::runtime_error(
                    "TabularDFT expects integer-degree phi1/phi2 samples"
                );
            }
            const T bindingEnergy =
                cfg.inputEnergyIsBinding
                ? inputEnergy
                : inputEnergy - T{2} * cfg.isolatedHelixEnergyEV;
            m_rawCurves[{
                positiveModulo(p1, 360),
                positiveModulo(p2, 360)
            }].push_back({
                R / cfg.L0,
                bindingEnergy / cfg.E0
            });
            ++accepted;
        }
        if (accepted == 0 || m_rawCurves.empty()) {
            throw std::runtime_error(
                "No rows found at requested zeta in " + filename
            );
        }
        for (auto& [key, curve] : m_rawCurves) {
            (void)key;
            std::sort(
                curve.begin(), curve.end(),
                [](const auto& a, const auto& b) { return a.r < b.r; }
            );
            if (curve.size() < 2) {
                throw std::runtime_error(
                    "Angular node has fewer than two R samples"
                );
            }
        }
    }
    // --------------------------------------------------------
    // Screw expansion:
    //
    // phi1' = phi1 + delta
    // phi2' = phi2 + screwDirection * delta
    // delta = 0, step, 2*step, ... < 360
    //
    // Multiple raw nodes mapping to the same expanded node are kept and
    // averaged later, matching the Python groupby(...).mean() behavior.
    // --------------------------------------------------------
    void buildScrewCoverage() {
        const int nShifts = 360 / m_screwStepDeg;
        for (const auto& [rawKey, curve] : m_rawCurves) {
            (void)curve;
            for (int k = 0; k < nShifts; ++k) {
                const int delta = k * m_screwStepDeg;
                const int p1 =
                    positiveModulo(rawKey.phi1 + delta, 360);
                const int p2 =
                    positiveModulo(
                        rawKey.phi2 + m_screwDirection * delta,
                        360
                    );
                auto& contributors =
                    m_expandedToRaw[angularIndex(p1, p2)];
                if (
                    std::find(
                        contributors.begin(),
                        contributors.end(),
                        rawKey
                    ) == contributors.end()
                ) {
                    contributors.push_back(rawKey);
                }
            }
        }
        m_missingAngularNodes = static_cast<std::size_t>(
            std::count_if(
                m_expandedToRaw.begin(),
                m_expandedToRaw.end(),
                [](const auto& x) { return x.empty(); }
            )
        );
    }
    // --------------------------------------------------------
    // Radial interpolation at one raw angular node.
    // Returns {U, dU/dr}.
    //
    // Below the first R:
    //   exponential repulsive continuation if the first slope is negative.
    //
    // Inside the table:
    //   linear interpolation.
    //
    // Beyond the final R:
    //   U = 0 and dU/dr = 0.
    // --------------------------------------------------------
    std::pair<T, T> evaluateRawCurve(
        const RadialCurve<T>& curve,
        T r
    ) const {
        if (r <= curve.front().r) {
            const auto& a = curve[0];
            const auto& b = curve[1];
            const T dr = b.r - a.r;
            if (std::abs(dr) <= std::numeric_limits<T>::epsilon()) {
                return {a.energy, T{0}};
            }
            const T slope = (b.energy - a.energy) / dr;
            if (slope < T{0}) {
                const T beta = T{12};
                const T x = a.r - r;
                const T amplitude = -slope / beta;
                const T expTerm = std::exp(beta * x);
                return {
                    a.energy + amplitude * (expTerm - T{1}),
                    slope * expTerm
                };
            }
            return {
                a.energy + slope * (r - a.r),
                slope
            };
        }
        if (r >= curve.back().r) {
            return {T{0}, T{0}};
        }
        const auto upper = std::upper_bound(
            curve.begin(), curve.end(), r,
            [](T value, const RadialPoint<T>& p) {
                return value < p.r;
            }
        );
        const std::size_t i1 = static_cast<std::size_t>(
            std::distance(curve.begin(), upper)
        );
        const auto& a = curve[i1 - 1];
        const auto& b = curve[i1];
        const T dr = b.r - a.r;
        if (std::abs(dr) <= std::numeric_limits<T>::epsilon()) {
            return {a.energy, T{0}};
        }
        const T slope = (b.energy - a.energy) / dr;
        return {
            a.energy + slope * (r - a.r),
            slope
        };
    }
    // --------------------------------------------------------
    // Evaluate one integer angular node after screw expansion.
    // If several raw nodes map here, average their radial curves.
    // --------------------------------------------------------
    std::pair<T, T> evaluateAngularNode(
        int phi1,
        int phi2,
        T r
    ) const {
        const auto& contributors =
            m_expandedToRaw[angularIndex(phi1, phi2)];
        if (contributors.empty()) {
            throw std::runtime_error(
                "TabularDFT reached an uncovered angular node"
            );
        }
        T energy = T{0};
        T dUdr = T{0};
        for (const auto& rawKey : contributors) {
            const auto [e, de] =
                evaluateRawCurve(m_rawCurves.at(rawKey), r);
            energy += e;
            dUdr += de;
        }
        const T norm =
            T{1} / static_cast<T>(contributors.size());
        return {energy * norm, dUdr * norm};
    }
};
// ============================================================
// Potential: owns EP, EA, OP and OA reference tables
// ============================================================
template<typename T>
class TabularDFTPotential {
public:
    explicit TabularDFTPotential(
        const std::string& referenceDirectory
    )
        : TabularDFTPotential(
            referenceDirectory,
            defaultTabularDFTConfig<T>()
        )
    {}

private:
    TabularDFTPotential(
        const std::string& referenceDirectory,
        const TabularDFTConfig<T>& cfg
    )
        : m_cutoff(cfg.cutoff),
          m_tables{
              makeTable(referenceDirectory, RefInteractionType::EP, cfg),
              makeTable(referenceDirectory, RefInteractionType::EA, cfg),
              makeTable(referenceDirectory, RefInteractionType::OP, cfg),
              makeTable(referenceDirectory, RefInteractionType::OA, cfg)
          }
    {
        std::cerr
            << "TabularDFT reduced units: E0=" << cfg.E0
            << " eV, L0=" << cfg.L0
            << " Angstrom, cutoff*=" << m_cutoff << '\n';
        std::cerr
            << "TabularDFT angular coverage missing nodes: "
            << "EP=" << table(RefInteractionType::EP).missingAngularNodes() << ", "
            << "EA=" << table(RefInteractionType::EA).missingAngularNodes() << ", "
            << "OP=" << table(RefInteractionType::OP).missingAngularNodes() << ", "
            << "OA=" << table(RefInteractionType::OA).missingAngularNodes() << '\n';
    }

public:
    T operator()(
        const ParticleOriented<T>& p1,
        const ParticleOriented<T>& p2,
        const Vec<T, 2>&,
        T r
    ) const {
        if (r <= T{0} || r > m_cutoff) return T{0};
        return table(getRefInteractionType(p1, p2))
            .evaluate(r, p1.phi, p2.phi)
            .energy;
    }
    std::pair<T, Vec<T, 3>> energyAndForces(
        const ParticleOriented<T>& p1,
        const ParticleOriented<T>& p2,
        const Vec<T, 2>& dr,
        T r
    ) const {
        if (r <= T{0} || r > m_cutoff) {
            return {T{0}, {{T{0}, T{0}, T{0}}}};
        }
        const auto value =
            table(getRefInteractionType(p1, p2))
                .evaluate(r, p1.phi, p2.phi);
        const T radialForce = -value.dUdr;
        const T fx = radialForce * dr[0] / r;
        const T fy = radialForce * dr[1] / r;
        // The force functor returns the generalized force acting on p2.
        const T torqueOn2 = -value.dUdPhi2;
        return {
            value.energy,
            {{fx, fy, torqueOn2}}
        };
    }
private:
    T m_cutoff;
    std::array<ScrewPeriodicReferenceTable<T>, 4> m_tables;
    const ScrewPeriodicReferenceTable<T>& table(
        RefInteractionType type
    ) const {
        return m_tables[static_cast<std::size_t>(type)];
    }
    static ScrewPeriodicReferenceTable<T> makeTable(
        const std::string& directory,
        RefInteractionType type,
        const TabularDFTConfig<T>& cfg
    ) {
        std::string filename = directory;
        if (!filename.empty() && filename.back() != '/') filename += '/';
        filename += interactionFilename(type);
        return ScrewPeriodicReferenceTable<T>(
            filename,
            cfg,
            inferScrewDirection(type)
        );
    }
};
// ============================================================
// Force wrapper
// ============================================================
template<typename TParticle, typename T = typename TParticle::value_type>
class TabularDFTForce {
public:
    using value_type = T;
    explicit TabularDFTForce(
        const std::string& referenceDirectory
    )
        : m_potential(referenceDirectory)
    {}
    Vec<T, 3> operator()(
        const TParticle& p1,
        const TParticle& p2,
        const Vec<T, 2>& dr,
        T r
    ) const {
        return m_potential.energyAndForces(p1, p2, dr, r).second;
    }
private:
    TabularDFTPotential<T> m_potential;
};
// ============================================================
// Factory for chiSurfMD
// ============================================================

template<typename Particle, typename SFINAE = void>
struct TabularDFT;

// ------------------------------------------------------------
template<typename T>
struct TabularDFT<ParticleOriented<T>> {
    static void initProgramOptions(
        po::options_description& desc
    ) {
        desc.add_options()
            (
                "dataReferenceDir",
                po::value<std::string>()->required(),
                "directory containing E_all_EP/EA/OP/OA.dat"
            )
            (
                "E0",
                po::value<T>()->default_value(T{0.4}),
                "reference energy in eV"
            )
            (
                "L0",
                po::value<T>()->default_value(T{7.0}),
                "reference length in Angstrom"
            );
    }

    // --------------------------------------------------------
    static auto force(
        const po::variables_map& vm
    ) {
        return TabularDFTForce<ParticleOriented<T>>(
            vm["dataReferenceDir"].as<std::string>()
        );
    }

    // --------------------------------------------------------
    static auto potential(
        const po::variables_map& vm
    ) {
        return TabularDFTPotential<T>(
            vm["dataReferenceDir"].as<std::string>()
        );
    }

    using ForceType =
        TabularDFTForce<ParticleOriented<T>>;
};
