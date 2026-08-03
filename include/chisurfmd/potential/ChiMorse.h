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

    T cutoff() const {
        return m_cutoff;
    }

    ChiMorseValue<T> evaluate(T r, T chi, T psi) const {
        if (r <= T{0} || r > m_cutoff) {
            return ChiMorseValue<T>{T{0}, T{0}, T{0}, T{0}};
        }

        const auto D = m_D.evaluate(chi, psi);
        const auto re = m_re.evaluate(chi, psi);
        const auto alpha = m_alpha.evaluate(chi, psi);

        const auto morse = evaluateMorse(
            r,
            D.value,
            re.value,
            alpha.value
        );

        ChiMorseValue<T> out;

        out.energy = morse.energy;
        out.dUdr = morse.dUdr;

        out.dUdChi =
            morse.dUdD * D.dChi
            + morse.dUdRe * re.dChi
            + morse.dUdAlpha * alpha.dChi;

        out.dUdPsi =
            morse.dUdD * D.dPsi
            + morse.dUdRe * re.dPsi
            + morse.dUdAlpha * alpha.dPsi;

        return out;
    }

private:
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
    const json& parameter
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

        FourierTerm<T> term;

        term.m = b.at("m").get<int>();
        term.n = b.at("n").get<int>();
        term.chiFunction = parseTrig(
            b.at("chi_function").get<std::string>()
        );
        term.psiFunction = parseTrig(
            b.at("psi_function").get<std::string>()
        );
        term.coefficient = coeffs[i];

        terms.push_back(term);
    }

    return FourierSurface2D<T>(std::move(terms));
}
// ------------------------------------------------------------
template<typename T>
AngularParameter<T> loadParameter(
    const json& basisTerms,
    const json& parameter
) {
    const std::string type =
        parameter.at("type").get<std::string>();

    if (type == "constant") {
        return AngularParameter<T>::constant(
            parameter.at("value").get<T>()
        );
    }

    if (type == "fourier") {
        return AngularParameter<T>::fourier(
            loadSurface<T>(basisTerms, parameter)
        );
    }

    throw std::runtime_error("Unknown parameter type: " + type);
}
// ------------------------------------------------------------
template<typename T>
ChiMorseModel<T> loadSingleChiMorseModel(
    const json& modelJson
) {
    const auto& basis =
        modelJson.at("basis_terms");

    const auto& params =
        modelJson.at("parameters");

    auto D =
        loadParameter<T>(
            basis,
            params.at("D")
        );

    auto re =
        loadParameter<T>(
            basis,
            params.at("re")
        );

    auto alpha =
        loadParameter<T>(
            basis,
            params.at("alpha")
        );

    const T cutoff =
        modelJson.at("cutoff").get<T>();

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
    const std::string& filename
) {
    std::ifstream file(filename);

    if (!file) {
        throw std::runtime_error(
            "Could not open chiMorse JSON file: "
            + filename
        );
    }

    json j;
    file >> j;

    const auto& interactions =
        j.at("interactions");

    const auto& factors =
        j.at("handedness_factor");

    return ChiMorseModels<T>{
        loadSingleChiMorseModel<T>(
            interactions.at("EP")
        ),
        loadSingleChiMorseModel<T>(
            interactions.at("EA")
        ),
        loadSingleChiMorseModel<T>(
            interactions.at("OP")
        ),
        loadSingleChiMorseModel<T>(
            interactions.at("OA")
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
        const auto type =
            getInteractionType(p1, p2);

        const auto& model =
            m_models.get(type);

        if (r <= T{0} || r > model.cutoff()) {
            return T{0};
        }

        const T h =
            m_models.getHandednessFactor(type);

        const T chi =
            p2.phi - h * p1.phi;

        const T psi =
            p1.phi + h * p2.phi;

        return model.evaluate(
            r,
            chi,
            psi
        ).energy;
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
        const auto type =
            getInteractionType(p1, p2);

        const auto& model =
            m_models.get(type);

        if (r <= T{0} || r > model.cutoff()) {
            return {{T{0}, T{0}, T{0}}};
        }

        const T h =
            m_models.getHandednessFactor(type);

        const T chi =
            p2.phi - h * p1.phi;

        const T psi =
            p1.phi + h * p2.phi;

        const auto value =
            model.evaluate(r, chi, psi);

        const T radialForce =
            -value.dUdr;

        const T fx =
            radialForce * dr[0] / r;

        const T fy =
            radialForce * dr[1] / r;

        const T torqueOn2 =
            -value.dUdChi
            -h * value.dUdPsi;

        return {{fx, fy, torqueOn2}};
    }

private:
    ChiMorseModels<T> m_models;
};

// ============================================================
// Factory struct for chiSurfMD
// ============================================================
namespace po = boost::program_options;
// ------------------------------------------------------------
template<typename Particle, typename SFINAE = void>
struct ChiMorse;
// ------------------------------------------------------------
template<typename T>
struct ChiMorse<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description& desc) {
        desc.add_options()
            (
                "chiMorseModel",
                po::value<std::string>()->required(),
                "path to chiMorse JSON model"
            );
    }

    static auto force(const po::variables_map& vm) {
        try {
            const auto filename =
                vm["chiMorseModel"].as<std::string>();

            auto models =
                loadChiMorseModels<T>(filename);

            return ChiMorseForce<ParticleOriented<T>>(
                std::move(models)
            );

        } catch (const boost::bad_any_cast& e) {
            std::cerr
                << "Error initializing ChiMorseForce: "
                << e.what()
                << std::endl;
            throw;
        }
    }

    static auto potential(const po::variables_map& vm) {
        try {
            const auto filename =
                vm["chiMorseModel"].as<std::string>();

            auto models =
                loadChiMorseModels<T>(filename);

            return ChiMorsePotential<ParticleOriented<T>>(
                std::move(models)
            );

        } catch (const boost::bad_any_cast& e) {
            std::cerr
                << "Error initializing ChiMorsePotential: "
                << e.what()
                << std::endl;
            throw;
        }
    }

    using ForceType = ChiMorseForce<ParticleOriented<T>>;
};
