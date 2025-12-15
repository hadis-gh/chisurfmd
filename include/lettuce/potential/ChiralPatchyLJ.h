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

// ================  Potential for ChiralPatchyLJ ================

template<typename TParticle, typename T = typename TParticle::value_type>
class ChiralPatchyLJPotential {
public:
    ChiralPatchyLJPotential(T epsilon_same, T epsilon_opp, 
                           T sigma_same, T sigma_opp, 
                           T cutoff, T patchSigma, int patchNums,
                           T chiral_angle_offset = 0.0)
        : m_epsilon_same(epsilon_same), m_epsilon_opp(epsilon_opp),
          m_sigma_same(sigma_same), m_sigma_opp(sigma_opp),
          m_cutoff(cutoff), m_patchSigma(patchSigma), 
          m_patchNums(patchNums), m_chiral_angle_offset(chiral_angle_offset)
    {
        m_twoPatchSigmaSq = 2.0 * patchSigma * patchSigma;
        m_patchAngles.reserve(m_patchNums);
        for (int a = 0; a < m_patchNums; a++) {
            m_patchAngles.push_back(2 * M_PI * a / m_patchNums);
        }
    }

    T operator()(const TParticle& p1, const TParticle& p2, 
                 const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return 0;

        // Determine if particles have same or opposite chirality
        T handedness_product = p1.h * p2.h;
        T epsilon, sigma;
        
        if (handedness_product > 0) { // Same chirality
            epsilon = m_epsilon_same;
            sigma = m_sigma_same;
        } else { // Opposite chirality
            epsilon = m_epsilon_opp;
            sigma = m_sigma_opp;
        }

        T baseLJ = lennardJones(R, sigma, epsilon);

        if (R < sigma) {
            return baseLJ;  // Pure repulsion in core region
        }

        const T dx = dr[0];
        const T dy = dr[1];
        const T phi_i = p1.phi;
        const T phi_j = p2.phi;

        T gamma_ij = std::atan2(dy, dx);  // Angle from i to j
        T gamma_ji = wrapAngle(gamma_ij + M_PI);  // Angle from j to i

        // Add chiral rotation to patch angles
        T chiral_rot_i = p1.h * m_chiral_angle_offset;
        T chiral_rot_j = p2.h * m_chiral_angle_offset;

        T sum_ei = 0;
        for (int l = 0; l < m_patchNums; ++l) {
            // FIXED: Add chiral rotation here
            T patchAng_i = m_patchAngles[l] + chiral_rot_i;
            T theta_i = wrapAngle(phi_i + patchAng_i - gamma_ij);
            sum_ei += std::exp(-theta_i * theta_i / m_twoPatchSigmaSq);
        }

        T sum_ej = 0;
        for (int m = 0; m < m_patchNums; ++m) {
            // FIXED: Add chiral rotation here
            T patch_ang_j = m_patchAngles[m] + chiral_rot_j;
            T theta_j = wrapAngle(phi_j + patch_ang_j - gamma_ji);
            sum_ej += std::exp(-theta_j * theta_j / m_twoPatchSigmaSq);
        }

        // Angular modulation factor
        T A = sum_ei * sum_ej;
        
        // Total potential: V = V_LJ * A
        return baseLJ * A;
    }

private:
    T m_epsilon_same, m_epsilon_opp, m_sigma_same, m_sigma_opp;
    T m_cutoff, m_patchSigma, m_chiral_angle_offset;
    T m_twoPatchSigmaSq;
    int m_patchNums;
    std::vector<T> m_patchAngles;
};

// ================  Force for ChiralPatchyLJ ================

template<typename TParticle, typename T = typename TParticle::value_type>
class ChiralPatchyLJForce {
public:
    using value_type = T;

    ChiralPatchyLJForce(T epsilon_same, T epsilon_opp, 
                       T sigma_same, T sigma_opp,
                       T cutoff, T patchSigma, int patchNums,
                       T chiral_angle_offset = 0.0)
        : m_epsilon_same(epsilon_same), m_epsilon_opp(epsilon_opp),
          m_sigma_same(sigma_same), m_sigma_opp(sigma_opp),
          m_cutoff(cutoff), m_patchSigma(patchSigma), 
          m_patchNums(patchNums), m_chiral_angle_offset(chiral_angle_offset)
    {
        m_twoPatchSigmaSq = 2.0 * patchSigma * patchSigma;
        m_patchAngles.reserve(m_patchNums);
        for (int a = 0; a < m_patchNums; a++) {
            m_patchAngles.push_back(2 * M_PI * a / m_patchNums);
        }
    }

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, 
                        const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};

        const T dx = dr[0];
        const T dy = dr[1];
        const T r2 = dx * dx + dy * dy;

        // Determine chiral parameters
        T handedness_product = p1.h * p2.h;
        T epsilon, sigma;
        
        if (handedness_product > 0) {  // Same chirality
            epsilon = m_epsilon_same;
            sigma = m_sigma_same;
        } else {  // Opposite chirality
            epsilon = m_epsilon_opp;
            sigma = m_sigma_opp;
        }

        T baseLJ = lennardJones(R, sigma, epsilon);
        T dVdR = lennardJonesDerivative(R, sigma, epsilon);

        if (R < sigma) {
            T fx = -dVdR * dx / R;
            T fy = -dVdR * dy / R;
            return {{fx, fy, 0}};  // No torque in core region
        }

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;

        T gamma_ij = std::atan2(dy, dx);
        T gamma_ji = wrapAngle(gamma_ij + M_PI);

        // Chirality-dependent patch rotations
        T chiral_rot_i = p1.h * m_chiral_angle_offset;
        T chiral_rot_j = p2.h * m_chiral_angle_offset;

        // Compute angular factors and their derivatives
        T sum_ei = 0, sum_ej = 0;
        T dsum_ei_dtheta = 0, dsum_ej_dtheta = 0;

        for (int l = 0; l < m_patchNums; ++l) {
            // FIXED: Add chiral rotation here
            T patchAng_i = m_patchAngles[l] + chiral_rot_i;
            T theta_i = wrapAngle(phi_i + patchAng_i - gamma_ij);
            T exp_i = std::exp(-theta_i * theta_i / m_twoPatchSigmaSq);
            sum_ei += exp_i;
            dsum_ei_dtheta += (-2.0 * theta_i / m_twoPatchSigmaSq) * exp_i;
        }

        for (int m = 0; m < m_patchNums; ++m) {
            // FIXED: Add chiral rotation here
            T patch_ang_j = m_patchAngles[m] + chiral_rot_j;
            T theta_j = wrapAngle(phi_j + patch_ang_j - gamma_ji);
            T exp_j = std::exp(-theta_j * theta_j / m_twoPatchSigmaSq);
            sum_ej += exp_j;
            dsum_ej_dtheta += (-2.0 * theta_j / m_twoPatchSigmaSq) * exp_j;
        }

        // Angular modulation factor
        T A = sum_ei * sum_ej;
        
        // Derivatives for torque calculation
        T dA_dtheta_i = sum_ej * dsum_ei_dtheta;  // ∂A/∂θ_i
        T dA_dtheta_j = sum_ei * dsum_ej_dtheta;  // ∂A/∂θ_j

        // Derivatives of gamma with respect to coordinates
        T dgamma_dx = -dy / r2;
        T dgamma_dy = dx / r2;

        // Chain rule for dA/dx and dA/dy
        T dA_dx = -(dA_dtheta_i + dA_dtheta_j) * dgamma_dx;
        T dA_dy = -(dA_dtheta_i + dA_dtheta_j) * dgamma_dy;

        // Spatial derivatives of R
        T dRdx = dx / R;
        T dRdy = dy / R;

        // Total force components: F = -∇V
        // ∇V = (∂V/∂R)(dR/dx) * A + V_LJ * (∂A/∂x)
        T fx = -(dVdR * dRdx * A + baseLJ * dA_dx);
        T fy = -(dVdR * dRdy * A + baseLJ * dA_dy);

        // Torques: τ = -∂V/∂φ
        T torque_i = -baseLJ * dA_dtheta_i;  // Torque on particle i
        T torque_j = -baseLJ * dA_dtheta_j;  // Torque on particle j

        // Return force on particle j and torque on particle j
        return {{fx, fy, torque_j}};
    }

private:
    T m_epsilon_same, m_epsilon_opp, m_sigma_same, m_sigma_opp;
    T m_cutoff, m_patchSigma, m_chiral_angle_offset;
    T m_twoPatchSigmaSq;
    int m_patchNums;
    std::vector<T> m_patchAngles;
};

// ================  Factory Struct ================

template<typename Particle, typename SFINAE = void>
struct ChiralPatchyLJ;

template<typename T>
struct ChiralPatchyLJ<ParticleOriented<T>>
{
    static void initProgramOptions(po::options_description &desc) {
        IsotropicLJ<ParticleDot<T>>::initProgramOptions(desc);

        desc.add_options()
            ("patchSigma",   po::value<T>()->default_value(0.3), 
                                   "Angular width of patches (radians)")
            ("patchNums",          po::value<int>()->default_value(2), 
                                   "Number of patches (symmetric)")
            ("epsilonSame",       po::value<T>()->default_value(1.0), 
                                   "Attraction strength for same chirality")
            ("epsilonOpp",        po::value<T>()->default_value(1.5), 
                                   "Attraction strength for opposite chirality")
            ("sigmaSame",         po::value<T>()->default_value(1.0), 
                                   "Equilibrium distance for same chirality")
            ("sigmaOpp",          po::value<T>()->default_value(0.9), 
                                   "Equilibrium distance for opposite chirality")
            ("chiralOffset",      po::value<T>()->default_value(0.2), 
                                   "Chiral twist of patches (radians)")
        ;
    }

    static auto force(const po::variables_map &vm) {
        return ChiralPatchyLJForce<ParticleOriented<T>>(
            vm["epsilonSame"].as<T>(), 
            vm["epsilonOpp"].as<T>(),
            vm["sigmaSame"].as<T>(), 
            vm["sigmaOpp"].as<T>(),
            vm["LJcutoff"].as<T>(),
            vm["patchSigma"].as<T>(),
            vm["patchNums"].as<int>(),
            vm["chiralOffset"].as<T>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return ChiralPatchyLJPotential<ParticleOriented<T>>(
            vm["epsilonSame"].as<T>(), 
            vm["epsilonOpp"].as<T>(),
            vm["sigmaSame"].as<T>(), 
            vm["sigmaOpp"].as<T>(),
            vm["LJcutoff"].as<T>(),
            vm["patchSigma"].as<T>(),
            vm["patchNums"].as<int>(),
            vm["chiralOffset"].as<T>()
        );
    }

    using ForceType = ChiralPatchyLJForce<ParticleOriented<T>>;    
};