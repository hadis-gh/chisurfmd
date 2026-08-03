#pragma once

#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <boost/program_options.hpp>
#include "chisurfmd/core/Vec.h"
#include "chisurfmd/core/ParticleDot.h"
#include "chisurfmd/core/ParticleOriented.h"
#include "chisurfmd/potential/IsotropicLJ.h"

namespace po = boost::program_options;

template<typename TParticle, typename T = typename TParticle::value_type>
class GeometricLJPotential {
public:
    GeometricLJPotential(T epsilon, T sigma, T cutoff, T patch_radius, int num_patches, T factor_patchy_lj)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), 
          m_patch_radius(patch_radius), m_num_patches(num_patches),
          m_factor_patchy_lj(factor_patchy_lj) {
        m_cutoff_sq = cutoff * cutoff;
    }

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return 0;

        T total_potential = lennardJones(R, m_sigma, m_epsilon);
        // T total_potential = 0;

        const T dx = dr[0];
        const T dy = dr[1];
        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T rho = m_patch_radius;

        // Patch-patch interaction (one-to-one mapping, n-th patch to n-th patch)
        for (int n = 0; n < m_num_patches; ++n) {
            const T theta_i_n = phi_i + 2.0 * M_PI * n / m_num_patches;
            const T theta_j_n = phi_j + 2.0 * M_PI * n / m_num_patches;

            const T delta_cos = rho * (std::cos(theta_j_n) - std::cos(theta_i_n));
            const T delta_sin = rho * (std::sin(theta_j_n) - std::sin(theta_i_n));

            const T r_patch_sq = (dx + delta_cos) * (dx + delta_cos) + 
                                 (dy + delta_sin) * (dy + delta_sin);

            const T r_patch = std::sqrt(r_patch_sq);
            
            total_potential += m_factor_patchy_lj * lennardJones(r_patch, m_sigma, m_epsilon);
        }

        return total_potential;     //Normalize based on numbers of patches!
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_patch_radius;
    T m_cutoff_sq;
    int m_num_patches;
    T m_factor_patchy_lj;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class GeometricLJForce {
public:
    GeometricLJForce(T epsilon, T sigma, T cutoff, T patch_radius, int num_patches, T factor_patchy_lj)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), 
          m_patch_radius(patch_radius), m_num_patches(num_patches),
          m_factor_patchy_lj(factor_patchy_lj) {
        m_cutoff_sq = cutoff * cutoff;
    }

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};

        const T dx = dr[0];
        const T dy = dr[1];
        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T rho = m_patch_radius;

        const T dRdx = dx / R;
        const T dRdy = dy / R;
        const T dVdR = lennardJonesDerivative(R, m_sigma, m_epsilon);
        // const T dVdR = 0;

        T fx = - dVdR * dRdx;
        T fy = - dVdR * dRdy; 
        T torque = 0;

        for (int n = 0; n < m_num_patches; ++n) {
            const T theta_i_n = phi_i + 2.0 * M_PI * n / m_num_patches;
            const T theta_j_n = phi_j + 2.0 * M_PI * n / m_num_patches;

            const T delta_cos = rho * (std::cos(theta_j_n) - std::cos(theta_i_n));
            const T delta_sin = rho * (std::sin(theta_j_n) - std::sin(theta_i_n));

            const T r_patch_dx = dx + delta_cos;
            const T r_patch_dy = dy + delta_sin;

            const T r_patch_sq = r_patch_dx * r_patch_dx + r_patch_dy * r_patch_dy;
            const T r_patch = std::sqrt(r_patch_sq);

            const T dVdr_patch = m_factor_patchy_lj * lennardJonesDerivative(r_patch, m_sigma, m_epsilon);

            const T dr_patchdx = r_patch_dx / r_patch;
            const T dr_patchdy = r_patch_dy / r_patch;

            fx += - dVdr_patch * dr_patchdx;
            fy += - dVdr_patch * dr_patchdy;

            const T ddelta_cos_dphi_j = rho * std::sin(theta_j_n);
            const T ddelta_sin_dphi_j = -rho * std::cos(theta_j_n);

            const T dr_patchdphi_j = (- ddelta_sin_dphi_j * r_patch_dy - ddelta_cos_dphi_j * r_patch_dx) / r_patch;
            torque += -dVdr_patch * dr_patchdphi_j;
        }
 
        return {{fx, fy, torque}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_patch_radius;
    T m_cutoff_sq;
    int m_num_patches;
    T m_factor_patchy_lj;
};

// ================================== Factory Struct for GeometricLJ ==================================

template<typename Particle, typename SFINAE = void>
struct GeometricLJ;

template<typename T>
struct GeometricLJ<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        IsotropicLJ<ParticleDot<T>>::initProgramOptions(desc);
        desc.add_options()
            ("factorPatchyLJ",      po::value<T>()->default_value(0.5),     "Strength of chiral interaction")
            ("patchNum",            po::value<int>()->default_value(1),     "Number of patches for geometric LJ")
            ("patchRadius",         po::value<T>()->default_value(0.5),     "Exclusion radius for geometric LJ")
        ;
    }

    static auto force(const po::variables_map &vm) {
        return GeometricLJForce<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["patchRadius"].as<T>(),
            vm["patchNum"].as<int>(),
            vm["factorPatchyLJ"].as<T>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return GeometricLJPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["patchRadius"].as<T>(),
            vm["patchNum"].as<int>(),
            vm["factorPatchyLJ"].as<T>()
        );
    }

    using ForceType = GeometricLJForce<ParticleOriented<T>>;
};