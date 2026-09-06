#pragma once

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <boost/program_options.hpp>
#include "nlohmann/json.hpp"

#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/ParticleOriented.h"

using json = nlohmann::json;
namespace po = boost::program_options;

// ============================================================
// Interaction Selection
// ============================================================
enum class InteractionType {
    EP,
    EA,
    OP,
    OA
};
// ------------------------------------------------------------
template<typename Particle>
InteractionType getInteractionType(
    const Particle& p1,
    const Particle& p2
) {
    const bool equalHandedness =
        p1.handedness == p2.handedness;

    const bool parallel =
        p1.alignment == p2.alignment;

    if (equalHandedness) {
        return parallel
            ? InteractionType::EP
            : InteractionType::EA;
    }

    return parallel
        ? InteractionType::OP
        : InteractionType::OA;
}
// ------------------------------------------------------------
// Canonical pair ordering
// 
//   EA    : +alignment particle first
//   OP/OA : +handedness particle first
//   EP    : already exchange-symmetric in the fitted basis
//
// This assumes that any global handedness/alignment transformation
// needed to construct the four reference tables has already been
// absorbed into their angle convention.
// ------------------------------------------------------------
template<typename Particle>
struct CanonicalPair {
    const Particle* first;
    const Particle* second;
    InteractionType type;
    bool swapped;
};
// ------------------------------------------------------------
template<typename Particle>
CanonicalPair<Particle> canonicalizePair(
    const Particle& p1,
    const Particle& p2
) {
    const auto type = getInteractionType(p1, p2);

    bool swap = false;

    switch (type) {
        case InteractionType::EP:
            break;

        case InteractionType::EA:
            // Same handedness, opposite alignment: +1 first.
            swap = p1.alignment < p2.alignment;
            break;

        case InteractionType::OP:
        case InteractionType::OA:
            // Opposite handedness: +1 first.
            swap = p1.handedness < p2.handedness;
            break;
    }

    if (swap) {
        return CanonicalPair<Particle>{
            &p2, &p1, type, true
        };
    }

    return CanonicalPair<Particle>{
        &p1, &p2, type, false
    };
}

// ============================================================
// Fourier evaluator
// ============================================================
enum class TrigType {
    Cos,
    Sin
// ------------------------------------------------------------
};
inline TrigType parseTrig(const std::string& s) {
    if (s == "cos") return TrigType::Cos;
    if (s == "sin") return TrigType::Sin;

    throw std::runtime_error("Unknown trig type: " + s);
}
// ------------------------------------------------------------
template<typename T>
struct FourierTerm {
    int m;
    int n;

    TrigType chiFunction;
    TrigType psiFunction;

    T coefficient;
};
// ------------------------------------------------------------
template<typename T>
struct FourierValue {
    T value;
    T dChi;
    T dPsi;
};
// ------------------------------------------------------------
template<typename T>
class FourierSurface2D {
public:
    FourierSurface2D() = default;

    explicit FourierSurface2D(std::vector<FourierTerm<T>> terms)
        : m_terms(std::move(terms)) {}

    FourierValue<T> evaluate(T chi, T psi) const {
        FourierValue<T> out{T{0}, T{0}, T{0}};

        for (const auto& term : m_terms) {
            const T m = static_cast<T>(term.m);
            const T n = static_cast<T>(term.n);

            const T mchi = m * chi;
            const T npsi = n * psi;

            T chiVal;
            T dChiVal;

            if (term.chiFunction == TrigType::Cos) {
                chiVal = std::cos(mchi);
                dChiVal = -m * std::sin(mchi);
            } else {
                chiVal = std::sin(mchi);
                dChiVal = m * std::cos(mchi);
            }

            T psiVal;
            T dPsiVal;

            if (term.psiFunction == TrigType::Cos) {
                psiVal = std::cos(npsi);
                dPsiVal = -n * std::sin(npsi);
            } else {
                psiVal = std::sin(npsi);
                dPsiVal = n * std::cos(npsi);
            }

            out.value += term.coefficient * chiVal * psiVal;
            out.dChi += term.coefficient * dChiVal * psiVal;
            out.dPsi += term.coefficient * chiVal * dPsiVal;
        }

        return out;
    }

    const std::vector<FourierTerm<T>>& terms() const {
        return m_terms;
    }

private:
    std::vector<FourierTerm<T>> m_terms;
};

// ============================================================
// Constant or Fourier parameter
// ============================================================
template<typename T>
class AngularParameter {
public:
    static AngularParameter constant(T value) {
        AngularParameter p;
        p.m_isConstant = true;
        p.m_constant = value;
        return p;
    }

    static AngularParameter fourier(FourierSurface2D<T> surface) {
        AngularParameter p;
        p.m_isConstant = false;
        p.m_surface = std::move(surface);
        return p;
    }

    FourierValue<T> evaluate(T chi, T psi) const {
        if (m_isConstant) {
            return FourierValue<T>{m_constant, T{0}, T{0}};
        }
        return m_surface.evaluate(chi, psi);
    }

    bool isConstant() const { return m_isConstant; }
    T constantValue() const { return m_constant; }
    const FourierSurface2D<T>& surface() const { return m_surface; }

private:
    bool m_isConstant = true;
    T m_constant = T{0};
    FourierSurface2D<T> m_surface;
};

// ============================================================
// Morse kernel
// ============================================================
template<typename T>
struct MorseValue {
    T energy;
    T dUdr;
    T dUdD;
    T dUdRe;
    T dUdAlpha;
};
// ------------------------------------------------------------
template<typename T>
MorseValue<T> evaluateMorse(T r, T D, T re, T alpha) {
    const T q = r - re;
    const T z = std::exp(-alpha * q);

    MorseValue<T> out;

    out.energy = D * (z * z - T{2} * z);

    out.dUdr = T{2} * D * alpha * z * (T{1} - z);

    out.dUdD = z * z - T{2} * z;

    out.dUdRe = -out.dUdr;

    out.dUdAlpha = T{2} * D * q * z * (T{1} - z);

    return out;
}

// ============================================================
// Full anisotropic Morse model
// ============================================================
template<typename T>
struct ChiMorseValue {
    T energy;
    T dUdr;
    T dUdChi;
    T dUdPsi;
};
// ------------------------------------------------------------
template<typename T>
class ChiMorseModel {
public:
    ChiMorseModel(
        AngularParameter<T> D,
        AngularParameter<T> re,
        AngularParameter<T> alpha,
        T cutoff
    )
        : m_D(std::move(D)),
          m_re(std::move(re)),
          m_alpha(std::move(alpha)),
          m_cutoff(cutoff) {}

    T cutoff() const { return m_cutoff; }

    ChiMorseValue<T> evaluate(T r, T chi, T psi) const {
        if (r <= T{0} || r > m_cutoff) {
            return {T{0}, T{0}, T{0}, T{0}};
        }

        FourierValue<T> D, re, alpha;
        evaluateAngular<true>(chi, psi, D, re, alpha);

        const auto morse = evaluateMorse(
            r, D.value, re.value, alpha.value
        );

        return {
            morse.energy,
            morse.dUdr,
            morse.dUdD * D.dChi
                + morse.dUdRe * re.dChi
                + morse.dUdAlpha * alpha.dChi,
            morse.dUdD * D.dPsi
                + morse.dUdRe * re.dPsi
                + morse.dUdAlpha * alpha.dPsi
        };
    }

    // Energy-only path: avoids all angular derivatives.
    T energy(T r, T chi, T psi) const {
        if (r <= T{0} || r > m_cutoff) return T{0};

        FourierValue<T> D, re, alpha;
        evaluateAngular<false>(chi, psi, D, re, alpha);

        const T q = r - re.value;
        const T z = std::exp(-alpha.value * q);
        return D.value * (z * z - T{2} * z);
    }

private:
    template<bool Derivatives>
    void evaluateAngular(
        T chi, T psi,
        FourierValue<T>& D,
        FourierValue<T>& re,
        FourierValue<T>& alpha
    ) const {
        const auto init = [](const AngularParameter<T>& p) {
            return p.isConstant()
                ? FourierValue<T>{p.constantValue(), T{0}, T{0}}
                : FourierValue<T>{T{0}, T{0}, T{0}};
        };

        D = init(m_D);
        re = init(m_re);
        alpha = init(m_alpha);

        const FourierSurface2D<T>* basis = nullptr;
        if (!m_D.isConstant()) basis = &m_D.surface();
        else if (!m_re.isConstant()) basis = &m_re.surface();
        else if (!m_alpha.isConstant()) basis = &m_alpha.surface();
        else return;

        const auto& terms = basis->terms();

        for (std::size_t i = 0; i < terms.size(); ++i) {
            const auto& b = terms[i];
            const T m = static_cast<T>(b.m);
            const T n = static_cast<T>(b.n);
            const T mc = m * chi;
            const T np = n * psi;

            T chiVal, psiVal, dChiVal = T{0}, dPsiVal = T{0};

            if constexpr (Derivatives) {
                if (b.chiFunction == TrigType::Cos) {
                    chiVal = std::cos(mc);
                    dChiVal = -m * std::sin(mc);
                } else {
                    chiVal = std::sin(mc);
                    dChiVal = m * std::cos(mc);
                }

                if (b.psiFunction == TrigType::Cos) {
                    psiVal = std::cos(np);
                    dPsiVal = -n * std::sin(np);
                } else {
                    psiVal = std::sin(np);
                    dPsiVal = n * std::cos(np);
                }
            } else {
                chiVal = (b.chiFunction == TrigType::Cos)
                    ? std::cos(mc) : std::sin(mc);
                psiVal = (b.psiFunction == TrigType::Cos)
                    ? std::cos(np) : std::sin(np);
            }

            const auto add = [&](const AngularParameter<T>& p,
                                 FourierValue<T>& out) {
                if (p.isConstant()) return;

                const T c = p.surface().terms()[i].coefficient;
                out.value += c * chiVal * psiVal;

                if constexpr (Derivatives) {
                    out.dChi += c * dChiVal * psiVal;
                    out.dPsi += c * chiVal * dPsiVal;
                }
            };

            add(m_D, D);
            add(m_re, re);
            add(m_alpha, alpha);
        }
    }

    AngularParameter<T> m_D;
    AngularParameter<T> m_re;
    AngularParameter<T> m_alpha;
    T m_cutoff;
};

template<typename T>
struct ChiMorseModels {
    ChiMorseModel<T> EP;
    ChiMorseModel<T> EA;
    ChiMorseModel<T> OP;
    ChiMorseModel<T> OA;

    T equalFactor;
    T oppositeFactor;

    const ChiMorseModel<T>& get(
        InteractionType type
    ) const {
        switch (type) {
            case InteractionType::EP: return EP;
            case InteractionType::EA: return EA;
            case InteractionType::OP: return OP;
            case InteractionType::OA: return OA;
        }

        throw std::runtime_error(
            "Unknown chiMorse interaction type"
        );
    }

    T getHandednessFactor(
        InteractionType type
    ) const {
        if (
            type == InteractionType::EP ||
            type == InteractionType::EA
        ) {
            return equalFactor;
        }

        return oppositeFactor;
    }
};

// ============================================================
// JSON loader
// ============================================================
template<typename T>
FourierSurface2D<T> loadSurface(
    const json& basisTerms,
    const json& parameter,
    T scale = T{1}
) {
    const auto coeffs =
        parameter.at("coefficients").template get<std::vector<T>>();

    if (coeffs.size() != basisTerms.size()) {
        throw std::runtime_error(
            "Coefficient size does not match basis size"
        );
    }

    std::vector<FourierTerm<T>> terms;
    terms.reserve(coeffs.size());

    for (std::size_t i = 0; i < coeffs.size(); ++i) {
        const auto& b = basisTerms.at(i);

        terms.push_back({
            b.at("m").get<int>(),
            b.at("n").get<int>(),
            parseTrig(b.at("chi_function").get<std::string>()),
            parseTrig(b.at("psi_function").get<std::string>()),
            coeffs[i] * scale
        });
    }

    return FourierSurface2D<T>(std::move(terms));
}
// ------------------------------------------------------------
template<typename T>
AngularParameter<T> loadParameter(
    const json& basisTerms,
    const json& parameter,
    T scale = T{1}
) {
    const auto type =
        parameter.at("type").get<std::string>();

    if (type == "constant") {
        return AngularParameter<T>::constant(
            parameter.at("value").get<T>() * scale
        );
    }

    if (type == "fourier") {
        return AngularParameter<T>::fourier(
            loadSurface<T>(basisTerms, parameter, scale)
        );
    }

    throw std::runtime_error(
        "Unknown parameter type: " + type
    );
}
// ------------------------------------------------------------
template<typename T>
ChiMorseModel<T> loadSingleChiMorseModel(
    const json& modelJson,
    T E0,
    T L0
) {
    const auto& basis =
        modelJson.at("basis_terms");

    const auto& params =
        modelJson.at("parameters");

    auto D = loadParameter<T>(
        basis,
        params.at("D"),
        T{1} / E0
    );

    auto re = loadParameter<T>(
        basis,
        params.at("re"),
        T{1} / L0
    );

    auto alpha = loadParameter<T>(
        basis,
        params.at("alpha"),
        L0
    );

    const T cutoff =
        modelJson.at("cutoff").get<T>() / L0;

    return ChiMorseModel<T>(
        std::move(D),
        std::move(re),
        std::move(alpha),
        cutoff
    );
}
// ------------------------------------------------------------
template<typename T>
ChiMorseModels<T> loadChiMorseModels(
    const std::string& filename,
    T E0,
    T L0
) {
    if (E0 <= T{0} || L0 <= T{0}) {
        throw std::runtime_error(
            "ChiMorse E0 and L0 must be positive"
        );
    }

    std::ifstream file(filename);

    if (!file) {
        throw std::runtime_error(
            "Could not open ChiMorse JSON file: " + filename
        );
    }

    json j;
    file >> j;

    const auto& interactions =
        j.at("interactions");

    const auto& factors =
        j.at("handedness_factor");

    return {
        loadSingleChiMorseModel<T>(
            interactions.at("EP"), E0, L0
        ),
        loadSingleChiMorseModel<T>(
            interactions.at("EA"), E0, L0
        ),
        loadSingleChiMorseModel<T>(
            interactions.at("OP"), E0, L0
        ),
        loadSingleChiMorseModel<T>(
            interactions.at("OA"), E0, L0
        ),
        factors.at("equal").get<T>(),
        factors.at("opposite").get<T>()
    };
}

// ============================================================
// Force and potential classes
// ============================================================
template<typename T>
struct ChiMorsePairForce {
    Vec<T, 2> forceOn1{{T{0}, T{0}}};

    T torqueOn1 = T{0};
    T torqueOn2 = T{0};
};
// ------------------------------------------------------------
template<
    typename TParticle,
    typename T = typename TParticle::value_type
>
class ChiMorsePotential {
public:
    explicit ChiMorsePotential(
        ChiMorseModels<T> models
    )
        : m_models(std::move(models)) {}

    T operator()(
        const TParticle& p1,
        const TParticle& p2,
        const Vec<T, 2>&,
        const T r
    ) const {
        const auto pair =
            canonicalizePair(p1, p2);

        const auto& model =
            m_models.get(pair.type);

        if (r <= T{0} || r > model.cutoff()) {
            return T{0};
        }

        const T h =
            m_models.getHandednessFactor(pair.type);

        const auto& c1 = *pair.first;
        const auto& c2 = *pair.second;

        const T chi =
            c2.phi - h * c1.phi;

        const T psi =
            c1.phi + h * c2.phi;

        return model.energy(r, chi, psi);
    }

private:
    ChiMorseModels<T> m_models;
};
// ------------------------------------------------------------
template<
    typename TParticle,
    typename T = typename TParticle::value_type
>
class ChiMorseForce {
public:
    using value_type = T;

    explicit ChiMorseForce(
        ChiMorseModels<T> models
    )
        : m_models(std::move(models)) {}

    Vec<T, 3> operator()(
        const TParticle& p1,
        const TParticle& p2,
        const Vec<T, 2>& dr,
        const T r
    ) const {
        const auto pair =
            canonicalizePair(p1, p2);

        const auto& model =
            m_models.get(pair.type);

        if (r <= T{0} || r > model.cutoff()) {
            return {{T{0}, T{0}, T{0}}};
        }

        const T h =
            m_models.getHandednessFactor(pair.type);

        const auto& c1 = *pair.first;
        const auto& c2 = *pair.second;

        const T chi =
            c2.phi - h * c1.phi;

        const T psi =
            c1.phi + h * c2.phi;

        const auto value =
            model.evaluate(r, chi, psi);

        const T radialForce =
            -value.dUdr;

        // dr always points from the actual p1 to the actual p2,
        // so these are the Cartesian force components on p2.
        const T fx =
            radialForce * dr[0] / r;

        const T fy =
            radialForce * dr[1] / r;

        // Torques of the canonical first and second particles.
        const T torqueOnCanonical1 =
            h * value.dUdChi
            - value.dUdPsi;

        const T torqueOnCanonical2 =
            -value.dUdChi
            -h * value.dUdPsi;

        // The force functor must return the torque on the actual p2.
        const T torqueOn2 = pair.swapped
            ? torqueOnCanonical1
            : torqueOnCanonical2;

        return {{fx, fy, torqueOn2}};
    }

private:
    ChiMorseModels<T> m_models;
};
// ============================================================
// Factory struct for ChiSurfMD
// ============================================================

template<typename Particle, typename SFINAE = void>
struct ChiMorse;


// ------------------------------------------------------------
template<typename T>
struct ChiMorse<ParticleOriented<T>>
{
    static void initProgramOptions(
        po::options_description& desc
    ) {
        desc.add_options()
            (
                "chiMorseModel",
                po::value<std::string>()->default_value("/home/hadis/chimorse/examples/models/chimorse_all.json"),
                "path to ChiMorse JSON model"
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


    static auto force(
        const po::variables_map& vm
    ) {
        auto models =
            loadChiMorseModels<T>(
                vm["chiMorseModel"].as<std::string>(),
                vm["E0"].as<T>(),
                vm["L0"].as<T>()
            );

        return ChiMorseForce<ParticleOriented<T>>(
            std::move(models)
        );
    }


    static auto potential(
        const po::variables_map& vm
    ) {
        auto models =
            loadChiMorseModels<T>(
                vm["chiMorseModel"].as<std::string>(),
                vm["E0"].as<T>(),
                vm["L0"].as<T>()
            );

        return ChiMorsePotential<ParticleOriented<T>>(
            std::move(models)
        );
    }


    using ForceType =
        ChiMorseForce<ParticleOriented<T>>;
};