#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"
#include <boost/program_options.hpp>
#include "IsotropicLJ.h"

namespace po = boost::program_options;

template<typename T>
inline T wrapAngle(const T& ang) {
    T result = std::fmod(ang + M_PI, 2.0 * M_PI);
    if (result < 0) result += 2.0 * M_PI;
    return result - M_PI;
}

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJPotential {
public:
    PatchyLJPotential(T epsilon, T sigma, T cutoff, T sigmaAngular, int patchNums)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), 
          m_sigmaAngular(sigmaAngular), m_patchNums(patchNums)
    {
        m_2sigmaAng_sq = 2.0 * sigmaAngular * sigmaAngular;
        m_patchAngles.reserve(m_patchNums);
        for (int a = 0; a < m_patchNums; a++) {
            m_patchAngles.push_back(2 * M_PI * a / m_patchNums);
        }
    }

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return 0;

        T baseLJ = lennardJones(R, m_sigma, m_epsilon);

        if (R < m_sigma) {
            return baseLJ;
        }

        const T dx = dr[0];
        const T dy = dr[1];
        const T phi_i = p1.phi;
        const T phi_j = p2.phi;

        T gamma_ij = std::atan2(dy, dx);
        T gamma_ji = wrapAngle(gamma_ij + M_PI);

        T sum_ei = 0;
        for (int l = 0; l < m_patchNums; ++l) {
            T patch_ang_i = m_patchAngles[l];
            T theta_i = wrapAngle(phi_i + patch_ang_i - gamma_ij);
            sum_ei += std::exp(-theta_i * theta_i / m_2sigmaAng_sq);
        }

        // Sum over patches for particle j
        T sum_ej = 0;
        for (int m = 0; m < m_patchNums; ++m) {
            T patch_ang_j = m_patchAngles[m];
            T theta_j = wrapAngle(phi_j + patch_ang_j - gamma_ji);
            sum_ej += std::exp(-theta_j * theta_j / m_2sigmaAng_sq);
        }

        // Angular modulation factor
        T A = sum_ei * sum_ej;

        // Total potential: V = V_LJ * A
        return baseLJ * A;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_2sigmaAng_sq;
    int m_patchNums;
    std::vector<T> m_patchAngles;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJForce {
public:
    using value_type = T;

    PatchyLJForce(T epsilon, T sigma, T cutoff, T sigmaAngular, int patchNums)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), 
          m_sigmaAngular(sigmaAngular), m_patchNums(patchNums)
    {
        m_2sigmaAng_sq = 2.0 * sigmaAngular * sigmaAngular;
        m_patchAngles.reserve(m_patchNums);
        for (int a = 0; a < m_patchNums; a++) {
            m_patchAngles.push_back(2 * M_PI * a / m_patchNums);
        }
    }

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};

        const T dx = dr[0];
        const T dy = dr[1];
        const T r2 = dx * dx + dy * dy;

        T baseLJ = lennardJones(R, m_sigma, m_epsilon);
        T dVdR = lennardJonesDerivative(R, m_sigma, m_epsilon);

        if (R < m_sigma) {
            T fx = -dVdR * dx / R;
            T fy = -dVdR * dy / R;
            return {{fx, fy, 0}};
        }

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;

        // Compute the angle from particle i to particle j
        T gamma_ij = std::atan2(dy, dx);
        // Angle from particle j to particle i
        T gamma_ji = wrapAngle(gamma_ij + M_PI);

        // Compute angular factors and their derivatives
        T sum_ei = 0, sum_ej = 0;
        T dsum_ei_dtheta = 0, dsum_ej_dtheta = 0;

        for (int l = 0; l < m_patchNums; ++l) {
            T patch_ang_i = m_patchAngles[l];
            T theta_i = wrapAngle(phi_i + patch_ang_i - gamma_ij);
            T exp_i = std::exp(-theta_i * theta_i / m_2sigmaAng_sq);
            sum_ei += exp_i;
            dsum_ei_dtheta += (-2.0 * theta_i / m_2sigmaAng_sq) * exp_i;
        }

        for (int m = 0; m < m_patchNums; ++m) {
            T patch_ang_j = m_patchAngles[m];
            T theta_j = wrapAngle(phi_j + patch_ang_j - gamma_ji);
            T exp_j = std::exp(-theta_j * theta_j / m_2sigmaAng_sq);
            sum_ej += exp_j;
            dsum_ej_dtheta += (-2.0 * theta_j / m_2sigmaAng_sq) * exp_j;
        }

        // Angular modulation factor and its partial derivatives
        T A = sum_ei * sum_ej;
        T dA_dtheta_i = sum_ej * dsum_ei_dtheta;  // ∂A/∂θ_i
        T dA_dtheta_j = sum_ei * dsum_ej_dtheta;  // ∂A/∂θ_j

        // Derivatives of gamma with respect to coordinates
        // ∂γ_ij/∂x = -dy/r², ∂γ_ij/∂y = dx/r²
        T dgamma_dx = -dy / r2;
        T dgamma_dy = dx / r2;

        // Chain rule for dA/dx and dA/dy
        // dA/dx = (∂A/∂θ_i)*(dθ_i/dγ)*(dγ/dx) + (∂A/∂θ_j)*(dθ_j/dγ)*(dγ/dx)
        // Since dθ_i/dγ = -1 and dθ_j/dγ = -1:
        T dA_dx = -(dA_dtheta_i + dA_dtheta_j) * dgamma_dx;
        T dA_dy = -(dA_dtheta_i + dA_dtheta_j) * dgamma_dy;

        // Spatial derivatives of R
        T dRdx = dx / R;
        T dRdy = dy / R;

        // Total force components: F = -∇V = -[ (dV/dR)*(dR/dx)*A + V*(dA/dx) ]
        T fx = -(dVdR * dRdx * A + baseLJ * dA_dx);
        T fy = -(dVdR * dRdy * A + baseLJ * dA_dy);

        // Torques: τ = -∂V/∂φ
        // ∂V/∂φ_i = V_LJ * (∂A/∂θ_i) * (∂θ_i/∂φ_i) = V_LJ * (∂A/∂θ_i) * 1
        T torque_i = -baseLJ * dA_dtheta_i;  // Torque on particle i
        T torque_j = -baseLJ * dA_dtheta_j;  // Torque on particle j

        // Return force on particle j and torque on particle j
        // Note: The force on particle i would be -fx, -fy with torque_i
        return {{fx, fy, torque_j}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigmaAngular;
    T m_2sigmaAng_sq;
    int m_patchNums;
    std::vector<T> m_patchAngles;
};

// ================  Factory Struct for PatchyLJ ================

template<typename Particle, typename SFINAE = void>
struct PatchyLJ;

template<typename T>
struct PatchyLJ<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        IsotropicLJ<ParticleDot<T>>::initProgramOptions(desc);

        desc.add_options()
            ("sigmaPatchyScale",   po::value<T>()->default_value(.262),              "anisotropic strength of potential")
            ("patchNums",          po::value<int>()->default_value(1),               "Number of patches for geometric LJ")
        ;
    }

    static auto force(const po::variables_map &vm) {
        return PatchyLJForce<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["sigmaPatchyScale"].as<T>(),
            vm["patchNums"].as<int>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return PatchyLJPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["sigmaPatchyScale"].as<T>(),
            vm["patchNums"].as<int>()
        );
    }

    static auto makePotential(const po::variables_map &vm) {
        return PatchyLJ<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["sigmaPatchyScale"].as<T>(),
            vm["patchNums"].as<int>()
        );
    }

    using ForceType = PatchyLJForce<ParticleOriented<T>>;
};