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
// 1. Interaction types and pair canonicalization
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

inline int inferScrewDirection(RefInteractionType type) {
    switch (type) {
        case RefInteractionType::EP:
        case RefInteractionType::EA: return +1;
        case RefInteractionType::OP:
        case RefInteractionType::OA: return -1;
    }
    throw std::runtime_error("Unknown TabularDFT interaction type");
}

template<typename T>
RefInteractionType getRefInteractionType(
    const ParticleOriented<T>& p1,
    const ParticleOriented<T>& p2
) {
    const bool equal    = p1.handedness == p2.handedness;
    const bool parallel = p1.alignment  == p2.alignment;
    if (equal) return parallel ? RefInteractionType::EP : RefInteractionType::EA;
    return parallel ? RefInteractionType::OP : RefInteractionType::OA;
}

template<typename T>
struct TabularCanonicalPair {
    const ParticleOriented<T>* first;
    const ParticleOriented<T>* second;
    RefInteractionType type;
    bool swapped;
};

template<typename T>
TabularCanonicalPair<T> canonicalizeTabularPair(
    const ParticleOriented<T>& p1,
    const ParticleOriented<T>& p2
) {
    const auto type = getRefInteractionType(p1, p2);
    bool swapped = false;
    switch (type) {
        case RefInteractionType::EP:
            break;
        case RefInteractionType::EA:
            swapped = p1.alignment  < p2.alignment;   // +1 first
            break;
        case RefInteractionType::OP:
        case RefInteractionType::OA:
            swapped = p1.handedness < p2.handedness;  // +1 first
            break;
    }
    return swapped
        ? TabularCanonicalPair<T>{&p2, &p1, type, true}
        : TabularCanonicalPair<T>{&p1, &p2, type, false};
}

// ============================================================
// 2. Configuration
// ============================================================

template<typename T>
struct TabularDFTConfig {
    T isolatedHelixEnergyEV;         // energy of one isolated helix (eV)
    bool inputEnergyIsBinding;       // true if files already store E_binding
    T E0;                            // energy scale (eV)
    T L0;                            // length scale (Angstrom)
    T cutoff;                        // requested cutoff in reduced units
    T zeta;                          // zeta slice to keep (Angstrom)
    T zetaTolerance;                 // accepted |zeta - zeta_target|
    int screwStepDeg;                // screw symmetry step (degrees)
    bool requireFullAngularCoverage; // abort if any 360x360 node is empty
    T shortRangeBeta;                // softness of below-r_min extension
    T switchStartFraction;           // fraction of cutoff where taper begins
};

namespace TabularDFTDefaults {
inline constexpr long double isolatedHelixEnergyEV = -6227.1749L;
inline constexpr bool        inputEnergyIsBinding  = false;
inline constexpr long double E0                    = 0.4L;
inline constexpr long double L0                    = 7.0L;
inline constexpr long double cutoffAngstrom        = 20.0L;
inline constexpr long double zeta                  = 0.0L;
inline constexpr long double zetaTolerance         = 1e-8L;
inline constexpr int         screwStepDeg          = 20;
inline constexpr bool        requireFullAngularCoverage = true;
inline constexpr long double shortRangeBeta        = 3.0L;
inline constexpr long double switchStartFraction   = 0.85L;
}

template<typename T>
TabularDFTConfig<T> defaultTabularDFTConfig() {
    return {
        static_cast<T>(TabularDFTDefaults::isolatedHelixEnergyEV),
        TabularDFTDefaults::inputEnergyIsBinding,
        static_cast<T>(TabularDFTDefaults::E0),
        static_cast<T>(TabularDFTDefaults::L0),
        static_cast<T>(TabularDFTDefaults::cutoffAngstrom /
                       TabularDFTDefaults::L0),
        static_cast<T>(TabularDFTDefaults::zeta),
        static_cast<T>(TabularDFTDefaults::zetaTolerance),
        TabularDFTDefaults::screwStepDeg,
        TabularDFTDefaults::requireFullAngularCoverage,
        static_cast<T>(TabularDFTDefaults::shortRangeBeta),
        static_cast<T>(TabularDFTDefaults::switchStartFraction)
    };
}

// ============================================================
// 3. Math helpers
// ============================================================

inline int positiveModulo(int x, int period) {
    const int r = x % period;
    return r < 0 ? r + period : r;
}

// Radians -> wrapped [0, 360) degrees.
template<typename T>
T radiansToWrappedDegrees(T angle) {
    constexpr T degPerRad = static_cast<T>(180.0L / 3.14159265358979323846L);
    angle = std::fmod(angle * degPerRad, T{360});
    return angle < T{0} ? angle + T{360} : angle;
}

// Quintic switch S(x), x in [0,1]:
//   S(0)=1, S'(0)=S''(0)=0     S(1)=0, S'(1)=S''(1)=0
// Tapers a function to zero without a force or force-derivative jump.
template<typename T>
inline T smoothSwitch(T x) {
    x = std::max(T{0}, std::min(T{1}, x));
    const T x3 = x * x * x;
    const T x4 = x3 * x;
    const T x5 = x4 * x;
    return T{1} - T{10} * x3 + T{15} * x4 - T{6} * x5;
}

template<typename T>
inline T smoothSwitchDerivative(T x) {
    x = std::max(T{0}, std::min(T{1}, x));
    const T x2 = x * x;
    const T x3 = x2 * x;
    const T x4 = x3 * x;
    return -T{30} * x2 + T{60} * x3 - T{30} * x4;
}

// Cubic Hermite on [0,1] with endpoint values (v0, v1) and endpoint
// tangents (m0, m1) per unit x. h = physical width of the segment.
// Returns { value, dvalue/dr } where r = r0 + x*h.
template<typename T>
struct Hermite1D { T value; T derivative; };

template<typename T>
Hermite1D<T> hermiteSegment(T x, T v0, T v1, T m0, T m1, T h) {
    const T x2 = x * x, x3 = x2 * x;
    const T h00 =  T{2} * x3 - T{3} * x2 + T{1};
    const T h10 =         x3 - T{2} * x2 + x;
    const T h01 = -T{2} * x3 + T{3} * x2;
    const T h11 =         x3 -        x2;
    const T dh00 =  T{6} * x2 - T{6} * x;
    const T dh10 =  T{3} * x2 - T{4} * x + T{1};
    const T dh01 = -T{6} * x2 + T{6} * x;
    const T dh11 =  T{3} * x2 - T{2} * x;
    Hermite1D<T> out;
    out.value      = h00 * v0 + h10 * h * m0 + h01 * v1 + h11 * h * m1;
    out.derivative = (dh00 * v0 + dh10 * h * m0 +
                      dh01 * v1 + dh11 * h * m1) / h;
    return out;
}

// Catmull-Rom: samples at x = -1, 0, 1, 2 (unit spacing); evaluate
// at a in [0,1] on the middle segment.
template<typename T>
Hermite1D<T> hermiteUniform(T a, T vm1, T v0, T v1, T v2) {
    const T m0 = T{0.5} * (v1 - vm1);
    const T m1 = T{0.5} * (v2 - v0);
    return hermiteSegment(a, v0, v1, m0, m1, T{1});
}

// ============================================================
// 4. Reference table (one per interaction type)
// ============================================================

template<typename T>
struct TabulatedValue {
    T energy  = T{0};
    T dUdr    = T{0};   // dU/dr*  (reduced)
    T dUdPhi1 = T{0};   // per radian
    T dUdPhi2 = T{0};   // per radian
};

template<typename T>
class ScrewPeriodicReferenceTable {
public:
    ScrewPeriodicReferenceTable(
        const std::string& filename,
        const TabularDFTConfig<T>& cfg,
        int screwDirection
    ) {
        validate(cfg);

        // Phase 1: read the raw file.
        RawTable raw = readRawTable(filename, cfg);

        // Phase 2: align all raw radial curves on a common r* grid.
        resampleToCommonGrid(raw);

        // Phase 3: screw-expand each raw node onto the 360x360 grid.
        buildExpandedGrid(raw, cfg.screwStepDeg, screwDirection);

        if (cfg.requireFullAngularCoverage && m_missingAngularNodes > 0) {
            throw std::runtime_error(
                "Incomplete screw-expanded angular grid in " + filename +
                ": missing " + std::to_string(m_missingAngularNodes) +
                " of 129600 nodes.");
        }

        // Effective cutoff is the smaller of the request and the data
        // extent; the smooth switch always ends at the effective cutoff.
        m_cutoff      = std::min(cfg.cutoff, m_radialGrid.back());
        m_switchStart = m_cutoff * cfg.switchStartFraction;
        m_beta        = cfg.shortRangeBeta;
    }

    std::size_t missingAngularNodes() const { return m_missingAngularNodes; }
    T cutoff()            const { return m_cutoff; }
    std::size_t nRadial() const { return m_radialGrid.size(); }
    T maxRadius()         const { return m_radialGrid.back(); }

    // --------------------------------------------------------
    // Evaluate at reduced r* and angles in radians.
    // Returns U, dU/dr*, dU/dphi1, dU/dphi2 (angular derivatives
    // are per radian; radial derivative is per reduced length).
    // --------------------------------------------------------
    TabulatedValue<T> evaluate(T r, T phi1Rad, T phi2Rad) const {
        TabulatedValue<T> out;
        if (r <= T{0} || r >= m_cutoff) return out;

        const T phi1Deg = radiansToWrappedDegrees(phi1Rad);
        const T phi2Deg = radiansToWrappedDegrees(phi2Rad);

        const int i1 = static_cast<int>(std::floor(phi1Deg));
        const int i2 = static_cast<int>(std::floor(phi2Deg));
        const T a1 = phi1Deg - static_cast<T>(i1);
        const T a2 = phi2Deg - static_cast<T>(i2);

        // 4x4 neighbourhood for Catmull-Rom.
        const int n1[4] = {
            positiveModulo(i1 - 1, 360), positiveModulo(i1,     360),
            positiveModulo(i1 + 1, 360), positiveModulo(i1 + 2, 360)
        };
        const int n2[4] = {
            positiveModulo(i2 - 1, 360), positiveModulo(i2,     360),
            positiveModulo(i2 + 1, 360), positiveModulo(i2 + 2, 360)
        };

        // Radial interpolation at each of the 16 nodes.
        T U_grid[4][4], dUdr_grid[4][4];
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 4; ++i) {
                const auto res = evaluateRadialAtNode(n1[i], n2[j], r);
                U_grid[i][j]    = res.first;
                dUdr_grid[i][j] = res.second;
            }
        }

        // Angular interpolation: rows (phi1), then columns (phi2).
        T U_row[4], dU1_row[4], dUdr_row[4];
        for (int j = 0; j < 4; ++j) {
            const auto hu = hermiteUniform(a1,
                U_grid[0][j], U_grid[1][j], U_grid[2][j], U_grid[3][j]);
            U_row[j]   = hu.value;
            dU1_row[j] = hu.derivative;
            const auto hd = hermiteUniform(a1,
                dUdr_grid[0][j], dUdr_grid[1][j],
                dUdr_grid[2][j], dUdr_grid[3][j]);
            dUdr_row[j] = hd.value;
        }
        const auto fu  = hermiteUniform(a2, U_row[0],
                                        U_row[1], U_row[2], U_row[3]);
        const auto fd1 = hermiteUniform(a2, dU1_row[0],
                                        dU1_row[1], dU1_row[2], dU1_row[3]);
        const auto fdr = hermiteUniform(a2, dUdr_row[0],
                                        dUdr_row[1], dUdr_row[2], dUdr_row[3]);

        constexpr T radPerDeg = static_cast<T>(3.14159265358979323846L / 180.0L);
        out.energy  = fu.value;
        out.dUdPhi1 = fd1.value / radPerDeg;
        out.dUdPhi2 = fu.derivative / radPerDeg;
        out.dUdr    = fdr.value;
        return out;
    }

private:
    // --------------------------------------------------------
    // Persistent state (all sizes fixed after construction).
    // --------------------------------------------------------
    std::vector<T> m_radialGrid;                // common r* samples
    std::vector<T> m_energy;                    // flat 360*360*n_r
    std::size_t    m_missingAngularNodes = 0;

    T m_beta        = T{0};
    T m_cutoff      = T{0};
    T m_switchStart = T{0};

    // --------------------------------------------------------
    // Construction scratch.
    // --------------------------------------------------------
    struct RawKey {
        int p1, p2;
        bool operator<(const RawKey& o) const {
            return p1 != o.p1 ? p1 < o.p1 : p2 < o.p2;
        }
    };
    using RawCurve = std::vector<std::pair<T, T>>;   // (r*, U*)
    using RawTable = std::map<RawKey, RawCurve>;

    static std::size_t angularIndex(int p1, int p2) {
        return static_cast<std::size_t>(p1 * 360 + p2);
    }

    static void validate(const TabularDFTConfig<T>& cfg) {
        if (!(cfg.E0 > T{0}))
            throw std::runtime_error("TabularDFT: E0 must be > 0");
        if (!(cfg.L0 > T{0}))
            throw std::runtime_error("TabularDFT: L0 must be > 0");
        if (cfg.screwStepDeg <= 0 || 360 % cfg.screwStepDeg != 0)
            throw std::runtime_error(
                "TabularDFT: screwStepDeg must divide 360");
    }

    // Read only the requested zeta slice, convert to binding energies
    // in reduced units, group by integer-degree (phi1, phi2).
    static RawTable readRawTable(
        const std::string& filename,
        const TabularDFTConfig<T>& cfg
    ) {
        std::ifstream in(filename);
        if (!in)
            throw std::runtime_error(
                "Could not open reference potential file: " + filename);

        RawTable raw;
        T phi1, phi2, zeta, R, e_in;
        std::size_t accepted = 0;
        while (in >> phi1 >> phi2 >> zeta >> R >> e_in) {
            if (std::abs(zeta - cfg.zeta) > cfg.zetaTolerance) continue;
            const int p1 = static_cast<int>(std::lround(phi1));
            const int p2 = static_cast<int>(std::lround(phi2));
            if (std::abs(phi1 - p1) > T{1e-6} ||
                std::abs(phi2 - p2) > T{1e-6})
                throw std::runtime_error(
                    "TabularDFT: non-integer phi samples in " + filename);
            const T binding = cfg.inputEnergyIsBinding
                ? e_in
                : e_in - T{2} * cfg.isolatedHelixEnergyEV;
            raw[{positiveModulo(p1, 360), positiveModulo(p2, 360)}]
                .emplace_back(R / cfg.L0, binding / cfg.E0);
            ++accepted;
        }
        if (accepted == 0 || raw.empty())
            throw std::runtime_error(
                "No rows at requested zeta in " + filename);
        for (auto& kv : raw) {
            auto& c = kv.second;
            std::sort(c.begin(), c.end(),
                [](const auto& a, const auto& b) { return a.first < b.first; });
            if (c.size() < 2)
                throw std::runtime_error(
                    "Angular node has fewer than two R samples");
        }
        return raw;
    }

    // Resample every raw curve onto a common r* grid (the densest
    // raw grid is used as reference; others are linearly interpolated).
    void resampleToCommonGrid(RawTable& raw) {
        std::size_t nMax = 0;
        const RawCurve* ref = nullptr;
        for (const auto& kv : raw)
            if (kv.second.size() > nMax) { nMax = kv.second.size(); ref = &kv.second; }

        m_radialGrid.resize(nMax);
        for (std::size_t i = 0; i < nMax; ++i)
            m_radialGrid[i] = (*ref)[i].first;

        for (auto& kv : raw) {
            const auto& src = kv.second;
            RawCurve dst(nMax);
            for (std::size_t i = 0; i < nMax; ++i) {
                const T q = m_radialGrid[i];
                dst[i].first = q;
                if (q <= src.front().first) { dst[i].second = src.front().second; continue; }
                if (q >= src.back().first)  { dst[i].second = src.back().second;  continue; }
                std::size_t lo = 0, hi = src.size() - 1;
                while (hi - lo > 1) {
                    const std::size_t mid = (lo + hi) / 2;
                    if (src[mid].first <= q) lo = mid; else hi = mid;
                }
                const T t = (q - src[lo].first) /
                            (src[hi].first - src[lo].first);
                dst[i].second = (T{1} - t) * src[lo].second +
                                        t  * src[hi].second;
            }
            kv.second = std::move(dst);
        }
    }

    // Screw expansion: raw node (phi1, phi2) contributes to
    //   (phi1 + delta, phi2 + s*delta)   for delta = k * step
    // for k = 0 .. 360/step-1. Duplicates are averaged.
    void buildExpandedGrid(
        const RawTable& raw,
        int screwStepDeg,
        int screwDirection
    ) {
        const std::size_t nR = m_radialGrid.size();
        const std::size_t N  = static_cast<std::size_t>(360) * 360 * nR;

        std::vector<T>   accum(N, T{0});
        std::vector<int> count(static_cast<std::size_t>(360) * 360, 0);

        const int nShifts = 360 / screwStepDeg;
        for (const auto& kv : raw) {
            const int rp1 = kv.first.p1;
            const int rp2 = kv.first.p2;
            const auto& curve = kv.second;
            for (int k = 0; k < nShifts; ++k) {
                const int delta = k * screwStepDeg;
                const int ep1   = positiveModulo(rp1 + delta, 360);
                const int ep2   = positiveModulo(rp2 + screwDirection * delta, 360);
                const std::size_t base = angularIndex(ep1, ep2) * nR;
                for (std::size_t i = 0; i < nR; ++i)
                    accum[base + i] += curve[i].second;
                ++count[angularIndex(ep1, ep2)];
            }
        }

        m_energy.assign(N, T{0});
        m_missingAngularNodes = 0;
        for (std::size_t a = 0; a < static_cast<std::size_t>(360) * 360; ++a) {
            if (count[a] == 0) { ++m_missingAngularNodes; continue; }
            const T inv = T{1} / static_cast<T>(count[a]);
            const std::size_t base = a * nR;
            for (std::size_t i = 0; i < nR; ++i)
                m_energy[base + i] = accum[base + i] * inv;
        }
    }

    // --------------------------------------------------------
    // Runtime: short-range continuation for r* < r*_min.
    // Slope-matched exponential with linear cap at depth 0.5 r*.
    // --------------------------------------------------------
    std::pair<T, T> shortRangeExtension(std::size_t base, T r) const {
        const T U0 = m_energy[base + 0];
        const T U1 = m_energy[base + 1];
        const T h0 = m_radialGrid[1] - m_radialGrid[0];
        const T slope = (U1 - U0) / h0;

        if (slope >= T{0})
            return {U0 + slope * (r - m_radialGrid[0]), slope};

        const T beta    = m_beta;
        const T x       = m_radialGrid[0] - r;   // > 0
        const T amp     = -slope / beta;
        const T capDepth = T{0.5};

        if (x > capDepth) {
            const T expCap = std::exp(beta * capDepth);
            const T Ucap   = U0 + amp * (expCap - T{1});
            const T dUcap  = slope * expCap;
            const T rCap   = m_radialGrid[0] - capDepth;
            return {Ucap + dUcap * (r - rCap), dUcap};
        }
        const T e = std::exp(beta * x);
        return {U0 + amp * (e - T{1}), slope * e};
    }

    // --------------------------------------------------------
    // Runtime: cubic Hermite in r* at one integer angular node,
    // with the smooth switch applied near the cutoff.
    // --------------------------------------------------------
    std::pair<T, T> evaluateRadialAtNode(int i1, int i2, T r) const {
        const std::size_t nR   = m_radialGrid.size();
        const std::size_t base = angularIndex(i1, i2) * nR;

        if (r <= m_radialGrid.front()) return shortRangeExtension(base, r);
        if (r >= m_radialGrid.back())  return {T{0}, T{0}};

        // Locate the segment [lo, lo+1] containing r.
        std::size_t lo = 0, hi = nR - 1;
        while (hi - lo > 1) {
            const std::size_t mid = (lo + hi) / 2;
            if (m_radialGrid[mid] <= r) lo = mid; else hi = mid;
        }
        const T r0 = m_radialGrid[lo];
        const T r1 = m_radialGrid[lo + 1];
        const T h  = r1 - r0;
        const T a  = (r - r0) / h;

        const T v0 = m_energy[base + lo];
        const T v1 = m_energy[base + lo + 1];

        auto u_at = [&](std::size_t i) -> T { return m_energy[base + i]; };

        // Catmull-Rom tangents from neighbouring samples.
        T m0, m1;
        if (lo == 0) {
            m0 = (v1 - v0) / h;
        } else {
            m0 = (u_at(lo + 1) - u_at(lo - 1)) /
                 (m_radialGrid[lo + 1] - m_radialGrid[lo - 1]);
        }
        if (lo + 2 >= nR) {
            m1 = (v1 - v0) / h;
        } else {
            m1 = (u_at(lo + 2) - u_at(lo)) /
                 (m_radialGrid[lo + 2] - m_radialGrid[lo]);
        }

        const auto hv = hermiteSegment(a, v0, v1, m0, m1, h);
        T U    = hv.value;
        T dUdr = hv.derivative;

        // C^2 taper to zero across [m_switchStart, m_cutoff].
        if (r > m_switchStart) {
            const T x  = (r - m_switchStart) / (m_cutoff - m_switchStart);
            const T S  = smoothSwitch(x);
            const T dS = smoothSwitchDerivative(x) / (m_cutoff - m_switchStart);
            dUdr = dUdr * S + U * dS;
            U    = U * S;
        }
        return {U, dUdr};
    }
};

// ============================================================
// 5. Potential: owns EP, EA, OP, OA tables
// ============================================================

template<typename T>
class TabularDFTPotential {
public:
    explicit TabularDFTPotential(const std::string& referenceDirectory)
        : TabularDFTPotential(referenceDirectory,
                              defaultTabularDFTConfig<T>()) {}

    T operator()(
        const ParticleOriented<T>& p1,
        const ParticleOriented<T>& p2,
        const Vec<T, 2>&,
        T r
    ) const {
        if (r <= T{0} || r >= m_cutoff) return T{0};
        const auto pair = canonicalizeTabularPair(p1, p2);
        return table(pair.type)
            .evaluate(r, pair.first->phi, pair.second->phi).energy;
    }

    std::pair<T, Vec<T, 3>> energyAndForces(
        const ParticleOriented<T>& p1,
        const ParticleOriented<T>& p2,
        const Vec<T, 2>& dr,
        T r
    ) const {
        if (r <= T{0} || r >= m_cutoff)
            return {T{0}, {{T{0}, T{0}, T{0}}}};

        const auto pair = canonicalizeTabularPair(p1, p2);
        const auto v = table(pair.type)
            .evaluate(r, pair.first->phi, pair.second->phi);

        const T F_r = -v.dUdr;
        const T fx  = F_r * dr[0] / r;
        const T fy  = F_r * dr[1] / r;

        // Force functor returns generalized force on caller's p2.
        const T torqueOn2 = pair.swapped ? -v.dUdPhi1 : -v.dUdPhi2;
        return {v.energy, {{fx, fy, torqueOn2}}};
    }

private:
    TabularDFTPotential(
        const std::string& dir,
        const TabularDFTConfig<T>& cfg
    )
        : m_cutoff(cfg.cutoff)
        , m_tables{
            makeTable(dir, RefInteractionType::EP, cfg),
            makeTable(dir, RefInteractionType::EA, cfg),
            makeTable(dir, RefInteractionType::OP, cfg),
            makeTable(dir, RefInteractionType::OA, cfg)
          }
    {
        std::cerr
            << "[TabularDFT] E0=" << cfg.E0
            << " eV, L0=" << cfg.L0
            << " Angstrom, cutoff*=" << cfg.cutoff
            << ", beta=" << cfg.shortRangeBeta
            << ", switch from " << cfg.switchStartFraction
            << " of cutoff\n";
        for (std::size_t i = 0; i < 4; ++i) {
            const auto t = static_cast<RefInteractionType>(i);
            const auto& tab = table(t);
            std::cerr
                << "[TabularDFT] " << interactionName(t)
                << ": n_r=" << tab.nRadial()
                << " r_max*=" << tab.maxRadius()
                << " missing_nodes=" << tab.missingAngularNodes()
                << '\n';
        }
    }

    T m_cutoff;
    std::array<ScrewPeriodicReferenceTable<T>, 4> m_tables;

    const ScrewPeriodicReferenceTable<T>& table(RefInteractionType t) const {
        return m_tables[static_cast<std::size_t>(t)];
    }

    static ScrewPeriodicReferenceTable<T> makeTable(
        const std::string& dir,
        RefInteractionType t,
        const TabularDFTConfig<T>& cfg
    ) {
        std::string fn = dir;
        if (!fn.empty() && fn.back() != '/') fn += '/';
        fn += interactionFilename(t);
        return ScrewPeriodicReferenceTable<T>(
            fn, cfg, inferScrewDirection(t));
    }
};

// ============================================================
// 6. Force wrapper
// ============================================================

template<typename TParticle, typename T = typename TParticle::value_type>
class TabularDFTForce {
public:
    using value_type = T;
    explicit TabularDFTForce(const std::string& dir) : m_potential(dir) {}

    Vec<T, 3> operator()(
        const TParticle& p1,
        const TParticle& p2,
        const Vec<T, 2>& dr,
        T r
    ) const {
        if (!std::isfinite(r))
            throw std::runtime_error(
                "TabularDFT: non-finite pair separation");
        const auto v = m_potential.energyAndForces(p1, p2, dr, r).second;
        if (!std::isfinite(v[0]) || !std::isfinite(v[1]) ||
            !std::isfinite(v[2]))
            throw std::runtime_error(
                "TabularDFT: non-finite force or torque");
        return v;
    }

private:
    TabularDFTPotential<T> m_potential;
};

// ============================================================
// 7. Program-options factory (CLI unchanged)
// ============================================================

template<typename Particle, typename SFINAE = void>
struct TabularDFT;

template<typename T>
struct TabularDFT<ParticleOriented<T>> {
    static void initProgramOptions(po::options_description& desc) {
        desc.add_options()
            ("dataReferenceDir",
             po::value<std::string>()->required(),
             "directory containing E_all_EP/EA/OP/OA.dat")
            ("E0", po::value<T>()->default_value(T{0.4}),
             "reference energy in eV")
            ("L0", po::value<T>()->default_value(T{7.0}),
             "reference length in Angstrom");
    }

    static auto force(const po::variables_map& vm) {
        return TabularDFTForce<ParticleOriented<T>>(
            vm["dataReferenceDir"].as<std::string>());
    }

    static auto potential(const po::variables_map& vm) {
        return TabularDFTPotential<T>(
            vm["dataReferenceDir"].as<std::string>());
    }

    using ForceType = TabularDFTForce<ParticleOriented<T>>;
};