#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "chisurfmd/potential/ChiMorse.h"
#include "chisurfmd/potential/TabularDFT.h"

namespace {
using Real = double;
using Clock = std::chrono::steady_clock;

constexpr Real pi = 3.14159265358979323846;
constexpr Real E0 = 0.4; // Same reduced units as mainMD.cpp / potential defaults
constexpr Real L0 = 7.0; // Angstrom

// EDIT THESE to a radial window well covered by your reference .dat files.
constexpr Real minR_A = 7.0;
constexpr Real maxR_A = 11.0;

// The benchmark uses the four interaction types in this order.
constexpr std::array<RefInteractionType, 4> types = {
    RefInteractionType::EP, RefInteractionType::EA,
    RefInteractionType::OP, RefInteractionType::OA
};
constexpr std::array<InteractionType, 4> chiTypes = {
    InteractionType::EP, InteractionType::EA,
    InteractionType::OP, InteractionType::OA
};

struct Configuration {
    std::size_t typeIndex; // EP, EA, OP, OA
    Real r;                // reduced distance, r (Angstrom) / L0
    Real phi1;             // radians
    Real phi2;             // radians
};

// Fixed random seed: both evaluators receive IDENTICAL configurations.
std::vector<Configuration> makeConfigurations(std::size_t n) {
    std::mt19937_64 rng(20261106);
    std::uniform_real_distribution<Real> radius(minR_A / L0, maxR_A / L0);
    std::uniform_real_distribution<Real> angle(0.0, 2.0 * pi);
    std::vector<Configuration> cases;
    cases.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        cases.push_back({i % 4, radius(rng), angle(rng), angle(rng)});
    }
    return cases;
}

template<class Evaluator>
void warmUp(const std::vector<Configuration>& cases, Evaluator&& evaluate) {
    volatile Real sink = 0.0;
    for (const auto& c : cases) sink += evaluate(c);
    (void)sink;
}

// Run full batches, then divide by the total number of evaluations.
// The checksum is used after timing so the compiler cannot discard the work.
template<class Evaluator>
Real timeNsPerCall(const std::vector<Configuration>& cases,
                   int passes, Evaluator&& evaluate, Real& checksum) {
    Real sum = 0.0;
    const auto start = Clock::now();
    for (int pass = 0; pass < passes; ++pass) {
        for (const auto& c : cases) sum += evaluate(c);
    }
    const auto finish = Clock::now();
    checksum += sum;
    const Real ns = std::chrono::duration<Real, std::nano>(finish - start).count();
    return ns / (static_cast<Real>(passes) * cases.size());
}

Real median(std::vector<Real> v) {
    std::sort(v.begin(), v.end());
    const auto m = v.size() / 2;
    return v.size() % 2 ? v[m] : (v[m - 1] + v[m]) / 2.0;
}

template<class TabEvaluator, class ChiEvaluator>
void compare(const std::string& operation,
             const std::vector<Configuration>& cases,
             int passes, int trials,
             TabEvaluator&& tabEval, ChiEvaluator&& chiEval,
             Real& checksum) {
    warmUp(cases, tabEval);
    warmUp(cases, chiEval);
    std::vector<Real> tabTimes, chiTimes;
    for (int trial = 0; trial < trials; ++trial) {
        if (trial % 2 == 0) {
            tabTimes.push_back(timeNsPerCall(cases, passes, tabEval, checksum));
            chiTimes.push_back(timeNsPerCall(cases, passes, chiEval, checksum));
        } else { // alternate to reduce systematic timing-order bias
            chiTimes.push_back(timeNsPerCall(cases, passes, chiEval, checksum));
            tabTimes.push_back(timeNsPerCall(cases, passes, tabEval, checksum));
        }
    }
    const Real tab = median(tabTimes);
    const Real chi = median(chiTimes);
    std::cout << std::left << std::setw(25) << operation
              << std::right << std::setw(15) << std::fixed << std::setprecision(1) << tab
              << std::setw(15) << chi
              << std::setw(16) << std::setprecision(2) << tab / chi << "\n";
}
}

int main(int argc, char** argv) {
    if (argc < 3 || argc > 6) {
        std::cerr << "Usage: " << argv[0]
                  << " REFERENCE_DATA_DIRECTORY CHIMORSE_MODEL.json"
                  << " [n_cases=2000] [passes=10] [trials=5]\n";
        return 1;
    }
    try {
        const std::string dataDir = argv[1];
        const std::string jsonPath = argv[2];
        const long long nInput = (argc > 3) ? std::stoll(argv[3]) : 2000;
        const int passes = (argc > 4) ? std::stoi(argv[4]) : 10;
        const int trials = (argc > 5) ? std::stoi(argv[5]) : 5;
        if (nInput < 4 || passes < 1 || trials < 1)
            throw std::runtime_error("n_cases must be >= 4; passes and trials must be positive");
        const auto n = static_cast<std::size_t>(nInput);

        // All file reading and model construction happens BEFORE timing.
        const auto cfg = defaultTabularDFTConfig<Real>();
        std::array<ScrewPeriodicReferenceTable<Real>, 4> tables = {
            ScrewPeriodicReferenceTable<Real>(dataDir + "/" + interactionFilename(types[0]), cfg, +1),
            ScrewPeriodicReferenceTable<Real>(dataDir + "/" + interactionFilename(types[1]), cfg, +1),
            ScrewPeriodicReferenceTable<Real>(dataDir + "/" + interactionFilename(types[2]), cfg, -1),
            ScrewPeriodicReferenceTable<Real>(dataDir + "/" + interactionFilename(types[3]), cfg, -1)
        };
        const auto models = loadChiMorseModels<Real>(jsonPath, E0, L0);
        for (const auto type : chiTypes) {
            if (maxR_A / L0 > models.get(type).cutoff())
                throw std::runtime_error("Selected radial window exceeds a ChiMorse cutoff");
        }
        if (maxR_A / L0 > cfg.cutoff)
            throw std::runtime_error("Selected radial window exceeds the tabulated cutoff");

        const auto cases = makeConfigurations(n);
        std::cout << "Radial window: " << minR_A << "-" << maxR_A
                  << " Angstrom | four interaction classes | " << n
                  << " configurations | " << passes << " passes | "
                  << trials << " trials\n";
        std::cout << "All times exclude loading and use the same configurations.\n\n";
        std::cout << std::left << std::setw(25) << "Operation"
                  << std::right << std::setw(15) << "Table ns/call"
                  << std::setw(15) << "Chi ns/call"
                  << std::setw(16) << "Table / Chi" << "\n";

        // 1. Energy + radial + angular derivatives: the FORCE/TORQUE kernel.
        // Table returns derivatives in phi1, phi2; ChiMorse in chi, psi.
        // Converting to the torque on particle 2 adds only a few operations.
        const auto tableDerivatives = [&](const Configuration& c) -> Real {
            const auto v = tables[c.typeIndex].evaluate(c.r, c.phi1, c.phi2);
            const Real torque2 = -v.dUdPhi2;
            return v.energy + 0.13 * v.dUdr + 0.013 * torque2;
        };
        const auto chiDerivatives = [&](const Configuration& c) -> Real {
            const auto type = chiTypes[c.typeIndex];
            const Real h = models.getHandednessFactor(type);
            const Real chi = c.phi2 - h * c.phi1;
            const Real psi = c.phi1 + h * c.phi2;
            const auto v = models.get(type).evaluate(c.r, chi, psi);
            const Real torque2 = -v.dUdChi - h * v.dUdPsi;
            return v.energy + 0.13 * v.dUdr + 0.013 * torque2;
        };

        // 2. ENERGY-ONLY: TabularDFT uses its existing evaluate() path,
        // while ChiMorse uses energy(). An optimizer may remove unused work.
        const auto tableEnergy = [&](const Configuration& c) -> Real {
            return tables[c.typeIndex].evaluate(c.r, c.phi1, c.phi2).energy;
        };
        const auto chiEnergy = [&](const Configuration& c) -> Real {
            const auto type = chiTypes[c.typeIndex];
            const Real h = models.getHandednessFactor(type);
            const Real chi = c.phi2 - h * c.phi1;
            const Real psi = c.phi1 + h * c.phi2;
            return models.get(type).energy(c.r, chi, psi);
        };

        Real checksum = 0.0;
        compare("Energy + derivatives", cases, passes, trials,
                tableDerivatives, chiDerivatives, checksum);
        compare("Energy only", cases, passes, trials,
                tableEnergy, chiEnergy, checksum);

        std::cout << "\nRatio > 1 means ChiMorse is faster for that operation.\n";
        std::cout << "Checksum (ignore value): " << std::setprecision(8)
                  << checksum << "\n";
        if (!std::isfinite(checksum))
            throw std::runtime_error("Non-finite result: check input data, units and radial window");
    } catch (const std::exception& e) {
        std::cerr << "Benchmark error: " << e.what() << "\n";
        return 2;
    }
    return 0;
}
