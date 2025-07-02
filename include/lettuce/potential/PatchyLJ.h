#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

template<typename T>
constexpr T wrapAngle(const T& ang) {
    return std::fmod(ang + M_PI, 2.0 * M_PI) - M_PI; 
}

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJPotential {
public:
    PatchyLJPotential(T epsilon, T sigma, T cutoff, T sigma_angular, unsigned int num_patches = 3, bool use_sum_mode = true)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), 
          m_sigma_angular(sigma_angular), m_num_patches(num_patches), m_use_sum_mode(use_sum_mode),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)),
          m_two_sigma_ang_sq(2.0 * sigma_angular * sigma_angular) {}

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return 0;

        // Base Lennard-Jones potential
        const T sr6 = m_sigma6 / std::pow(R, 6);
        const T sr12 = sr6 * sr6;
        const T LJ_pot = 4.0 * m_epsilon * (sr12 - sr6);

        // Pure repulsion region (r < sigma)
        if (R < m_sigma) {
            return LJ_pot;
        }

        // Calculate angular modulation factor
        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T gamma_ij = std::atan2(dr[1], dr[0]);
        const T gamma_ji = gamma_ij + M_PI;  // gamma_ji = gamma_ij + π

        T angular_factor = 0;

        if (m_use_sum_mode) {
            // Sum over all patch combinations
            for (int n = 0; n < m_num_patches; ++n) {
                const T psi_i_n = phi_i + 2.0 * M_PI * n / m_num_patches;
                const T psi_j_n = phi_j + 2.0 * M_PI * n / m_num_patches;
                
                const T theta_i_n = wrapAngle(psi_i_n - gamma_ij);
                const T theta_j_n = wrapAngle(psi_j_n - gamma_ji);
                
                const T exp_i = std::exp(-theta_i_n * theta_i_n / m_two_sigma_ang_sq);
                const T exp_j = std::exp(-theta_j_n * theta_j_n / m_two_sigma_ang_sq);
                
                angular_factor += exp_i * exp_j;
            }
        } else {
            // Maximum over all patch combinations (original patchy colloid model)
            T max_factor = 0;
            for (int n = 0; n < m_num_patches; ++n) {
                for (int m = 0; m < m_num_patches; ++m) {
                    const T psi_i_n = phi_i + 2.0 * M_PI * n / m_num_patches;
                    const T psi_j_m = phi_j + 2.0 * M_PI * m / m_num_patches;
                    
                    const T theta_i_n = wrapAngle(psi_i_n - gamma_ij);
                    const T theta_j_m = wrapAngle(psi_j_m - gamma_ji);
                    
                    const T exp_i = std::exp(-theta_i_n * theta_i_n / m_two_sigma_ang_sq);
                    const T exp_j = std::exp(-theta_j_m * theta_j_m / m_two_sigma_ang_sq);
                    
                    max_factor = std::max(max_factor, exp_i * exp_j);
                }
            }
            angular_factor = max_factor;
        }

        return LJ_pot * angular_factor;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigma_angular;
    unsigned int m_num_patches;
    bool m_use_sum_mode;
    T m_sigma6, m_sigma12, m_two_sigma_ang_sq;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJForce {
public:
    PatchyLJForce(T epsilon, T sigma, T cutoff, T sigma_angular, unsigned int num_patches = 3, bool use_sum_mode = true)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), 
          m_sigma_angular(sigma_angular), m_num_patches(num_patches), m_use_sum_mode(use_sum_mode),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)),
          m_two_sigma_ang_sq(2.0 * sigma_angular * sigma_angular) {}

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};

        // Base Lennard-Jones potential and force
        const T sr6 = m_sigma6 / std::pow(R, 6);
        const T sr12 = sr6 * sr6;
        const T LJ_pot = 4.0 * m_epsilon * (sr12 - sr6);
        const T dLJ_dR = 24.0 * m_epsilon / R * (2.0 * sr12 - sr6);

        // Pure repulsion region (r < sigma) - no angular dependence
        if (R < m_sigma) {
            const T force_mag = dLJ_dR;
            return {{force_mag * dr[0] / R, force_mag * dr[1] / R, 0}};
        }

        // Calculate angular terms and derivatives
        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T gamma_ij = std::atan2(dr[1], dr[0]);
        const T gamma_ji = gamma_ij + M_PI;
        
        const T R_sq = R * R;
        const T dgamma_dx = -dr[1] / R_sq;  // ∂γ/∂x
        const T dgamma_dy = dr[0] / R_sq;   // ∂γ/∂y

        T angular_factor = 0;
        T dAngular_dx = 0, dAngular_dy = 0, dAngular_dphi_i = 0;

        if (m_use_sum_mode) {
            // Sum mode - differentiate the sum
            for (int n = 0; n < m_num_patches; ++n) {
                const T psi_i_n = phi_i + 2.0 * M_PI * n / m_num_patches;
                const T psi_j_n = phi_j + 2.0 * M_PI * n / m_num_patches;
                
                const T theta_i_n = wrapAngle(psi_i_n - gamma_ij);
                const T theta_j_n = wrapAngle(psi_j_n - gamma_ji);
                
                const T exp_i = std::exp(-theta_i_n * theta_i_n / m_two_sigma_ang_sq);
                const T exp_j = std::exp(-theta_j_n * theta_j_n / m_two_sigma_ang_sq);
                const T term = exp_i * exp_j;
                
                angular_factor += term;
                
                // Derivatives
                const T dexp_i_dtheta_i = -theta_i_n * exp_i / m_sigma_angular / m_sigma_angular;
                const T dexp_j_dtheta_j = -theta_j_n * exp_j / m_sigma_angular / m_sigma_angular;
                
                // Chain rule: ∂Λ/∂x = ∂Λ/∂θᵢ * ∂θᵢ/∂γ * ∂γ/∂x + ∂Λ/∂θⱼ * ∂θⱼ/∂γ * ∂γ/∂x
                dAngular_dx += (dexp_i_dtheta_i * exp_j * dgamma_dx + 
                               exp_i * dexp_j_dtheta_j * (-dgamma_dx));  // Note: ∂γⱼᵢ/∂x = -∂γᵢⱼ/∂x
                
                dAngular_dy += (dexp_i_dtheta_i * exp_j * dgamma_dy + 
                               exp_i * dexp_j_dtheta_j * (-dgamma_dy));
                
                dAngular_dphi_i += dexp_i_dtheta_i * exp_j;  // ∂θᵢ/∂φᵢ = 1
            }
        } else {
            // Maximum mode - more complex derivative calculation
            T max_factor = 0;
            int best_n = 0, best_m = 0;
            
            // Find the maximum contributing patch pair
            for (int n = 0; n < m_num_patches; ++n) {
                for (int m = 0; m < m_num_patches; ++m) {
                    const T psi_i_n = phi_i + 2.0 * M_PI * n / m_num_patches;
                    const T psi_j_m = phi_j + 2.0 * M_PI * m / m_num_patches;
                    
                    const T theta_i_n = wrapAngle(psi_i_n - gamma_ij);
                    const T theta_j_m = wrapAngle(psi_j_m - gamma_ji);
                    
                    const T exp_i = std::exp(-theta_i_n * theta_i_n / m_two_sigma_ang_sq);
                    const T exp_j = std::exp(-theta_j_m * theta_j_m / m_two_sigma_ang_sq);
                    const T factor = exp_i * exp_j;
                    
                    if (factor > max_factor) {
                        max_factor = factor;
                        best_n = n;
                        best_m = m;
                    }
                }
            }
            
            angular_factor = max_factor;
            
            // Calculate derivatives for the best patch pair
            if (max_factor > 0) {
                const T psi_i_best = phi_i + 2.0 * M_PI * best_n / m_num_patches;
                const T psi_j_best = phi_j + 2.0 * M_PI * best_m / m_num_patches;
                
                const T theta_i_best = wrapAngle(psi_i_best - gamma_ij);
                const T theta_j_best = wrapAngle(psi_j_best - gamma_ji);
                
                const T exp_i = std::exp(-theta_i_best * theta_i_best / m_two_sigma_ang_sq);
                const T exp_j = std::exp(-theta_j_best * theta_j_best / m_two_sigma_ang_sq);
                
                const T dexp_i_dtheta_i = -theta_i_best * exp_i / m_sigma_angular / m_sigma_angular;
                const T dexp_j_dtheta_j = -theta_j_best * exp_j / m_sigma_angular / m_sigma_angular;
                
                dAngular_dx = dexp_i_dtheta_i * exp_j * dgamma_dx + exp_i * dexp_j_dtheta_j * (-dgamma_dx);
                dAngular_dy = dexp_i_dtheta_i * exp_j * dgamma_dy + exp_i * dexp_j_dtheta_j * (-dgamma_dy);
                dAngular_dphi_i = dexp_i_dtheta_i * exp_j;
            }
        }

        // Total force and torque
        const T force_x = -(dLJ_dR * angular_factor * dr[0] / R + LJ_pot * dAngular_dx);
        const T force_y = -(dLJ_dR * angular_factor * dr[1] / R + LJ_pot * dAngular_dy);
        const T torque = -LJ_pot * dAngular_dphi_i;

        return {{force_x, force_y, torque}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_sigma_angular;
    unsigned int m_num_patches;
    bool m_use_sum_mode;
    T m_sigma6, m_sigma12, m_two_sigma_ang_sq;
};