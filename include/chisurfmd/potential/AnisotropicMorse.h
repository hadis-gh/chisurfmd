#pragma once

#include <vector>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <boost/program_options.hpp>
#include "core/Vec.h"
#include "core/ParticleOriented.h"

namespace po = boost::program_options;


// =============================== Helper for Fourier Evaluation ===============================

template<typename T>
T evaluateFourierSeries(
    const T chi, const T psi,
    const std::vector<T>& coeffs,
    int h_chi, int h_psi,
    bool symm_chi, int k0) {

    std::vector<T> terms;
    terms.push_back(1.0); // Constant term

    // Pure chi terms
    for (int m = 1; m <= h_chi; ++m) {
        terms.push_back(std::cos(m * chi));
        if (!symm_chi) {
            terms.push_back(std::sin(m * chi));
        }
    }

    // Pure psi terms
    for (int j = 1; j <= h_psi; ++j) {
        const T n_psi = k0 * j * psi;
        terms.push_back(std::cos(n_psi));
        terms.push_back(std::sin(n_psi));
    }

    // Coupled chi-psi terms
    for (int m = 1; m <= h_chi; ++m) {
        T cos_m_chi = std::cos(m * chi);
        for (int j = 1; j <= h_psi; ++j) {
            const T n_psi = k0 * j * psi;
            T cos_n_psi = std::cos(n_psi);
            T sin_n_psi = std::sin(n_psi);

            terms.push_back(cos_m_chi * cos_n_psi);
            terms.push_back(cos_m_chi * sin_n_psi);
            if (!symm_chi) {
                T sin_m_chi = std::sin(m * chi);
                terms.push_back(sin_m_chi * cos_n_psi);
                terms.push_back(sin_m_chi * sin_n_psi);
            }
        }
    }

    if (coeffs.size() != terms.size()) {
        throw std::invalid_argument( "missmatch between coefficients size and terms numbers of Fourier expansion!" );
    }

    return std::inner_product(terms.begin(), terms.end(), coeffs.begin(), 0.0);
}


// =================================== Anisotropic Morse Potential ===================================

template<typename TParticle, typename T = typename TParticle::value_type>
class AnisotropicMorsePotential {
public:
    AnisotropicMorsePotential(T alpha, T cutoff, int h_chi, int h_psi, bool symm_chi, int screw_symm,
                              const std::vector<T>& D_coeffs, const std::vector<T>& re_coeffs)
        : m_alpha(alpha), m_cutoff(cutoff), m_h_chi(h_chi), m_h_psi(h_psi),
          m_symm_chi(symm_chi), m_D_coeffs(D_coeffs), m_re_coeffs(re_coeffs) {
        m_k0 = static_cast<int>(round(360.0 / (2.0 * screw_symm)));
        m_cutoff_sq = cutoff * cutoff;
    }

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return 0;

        const T chi = p1.phi - p2.phi;
        const T psi = std::atan2(dr[1], dr[0]);

        const T D = evaluateFourierSeries(chi, psi, m_D_coeffs, m_h_chi, m_h_psi, m_symm_chi, m_k0);
        const T re = evaluateFourierSeries(chi, psi, m_re_coeffs, m_h_chi, m_h_psi, m_symm_chi, m_k0);
        
        const T exponent_term = std::exp(-m_alpha * (R - re));
        return D * (exponent_term * exponent_term - 2 * exponent_term);
    }

private:
    T m_alpha, m_cutoff;
    int m_h_chi, m_h_psi, m_k0;
    bool m_symm_chi;
    std::vector<T> m_D_coeffs, m_re_coeffs;
    T m_cutoff_sq;
};

// =================================== Anisotropic Morse Force ===================================

template<typename TParticle, typename T = typename TParticle::value_type>
class AnisotropicMorseForce {
public:
    AnisotropicMorseForce(T alpha, T cutoff, int h_chi, int h_psi, bool symm_chi, int screw_symm,
                              const std::vector<T>& D_coeffs, const std::vector<T>& re_coeffs)
        : m_alpha(alpha), m_cutoff(cutoff), m_h_chi(h_chi), m_h_psi(h_psi),
          m_symm_chi(symm_chi), m_D_coeffs(D_coeffs), m_re_coeffs(re_coeffs) {
        m_k0 = static_cast<int>(round(360.0 / (2.0 * screw_symm)));
        m_cutoff_sq = cutoff * cutoff;
    }

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};
        
        // Place holder for the Force!

        // 1. A function `evaluateFourierSeries_Derivatives` that computes dD/dchi, dD/dpsi, etc.
        // 2. The chain rule to get dU/dx, dU/dy, and dU/dphi_i.
        //    dU/dx = (∂U/∂R)(∂R/∂x) + (∂U/∂psi)(∂psi/∂x)
        //    dU/dy = (∂U/∂R)(∂R/∂y) + (∂U/∂psi)(∂psi/∂y)
        //    dU/dphi_i = (∂U/∂chi)(∂chi/∂phi_i)

        return {{0, 0, 0}}; 
    }

private:
    T m_alpha, m_cutoff;
    int m_h_chi, m_h_psi, m_k0;
    bool m_symm_chi;
    std::vector<T> m_D_coeffs, m_re_coeffs;
    T m_cutoff_sq;
};


// =============================== Factory Struct for AnisotropicMorse ===============================
template<typename Particle, typename SFINAE = void>
struct AnisotropicMorse;

template<typename T>
struct AnisotropicMorse<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        desc.add_options()
            ("morseAlpha",      po::value<T>()->default_value(1.1),     "Alpha parameter for Morse potential")
            ("morseCutoff",     po::value<T>()->default_value(2.5),     "Cutoff distance for Morse potential")
            ("morseHChi",       po::value<int>()->default_value(6),     "Max harmonic m for chi")
            ("morseHPsi",       po::value<int>()->default_value(2),     "Max harmonic j for psi")
            ("morseSymmChi",    po::value<bool>()->default_value(true), "Symmetry in chi direction")
            ("morseScrewSymm",  po::value<int>()->default_value(20),    "Screw symmetry angle in degrees")
        ;
    }
    
    static auto force(const po::variables_map &vm) {
        // --- PLACEHOLDER COEFFICIENTS ---
        // The size must match what `evaluateFourierSeries` expects.
        std::vector<T> d_coeffs = {1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6}; 
        std::vector<T> re_coeffs = {1.1, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06}; 

        return AnisotropicMorseForce<ParticleOriented<T>>(
            vm["morseAlpha"].as<T>(),
            vm["morseCutoff"].as<T>(),
            vm["morseHChi"].as<int>(),
            vm["morseHPsi"].as<int>(),
            vm["morseSymmChi"].as<bool>(),
            vm["morseScrewSymm"].as<int>(),
            d_coeffs,
            re_coeffs
        );
    }

    static auto potential(const po::variables_map &vm) {

        std::vector<T> d_coeffs = {1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6}; 
        std::vector<T> re_coeffs = {1.1, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06}; 

        return AnisotropicMorsePotential<ParticleOriented<T>>(
            vm["morseAlpha"].as<T>(),
            vm["morseCutoff"].as<T>(),
            vm["morseHChi"].as<int>(),
            vm["morseHPsi"].as<int>(),
            vm["morseSymmChi"].as<bool>(),
            vm["morseScrewSymm"].as<int>(),
            d_coeffs,
            re_coeffs
        );
    }

    using ForceType = AnisotropicMorseForce<ParticleOriented<T>>;
};