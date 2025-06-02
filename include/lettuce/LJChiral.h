#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/ParticleDot.h"
#include "lettuce/ParticleOriented.h"

template<typename TParticle, typename T = typename TParticle::value_type>
class LJChiralForce {
public:
    LJChiralForce(T epsilon, T sigma, T cutoff, T particle_radius)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_particle_radius(particle_radius),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)) {}

    Vec<T, 2> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {
        if (r == 0 || r > m_cutoff) return {{0, 0}};

        T phi_i = p1.phi;
        T phi_j = p2.phi;
        T theta = std::atan2(dr[1], dr[0]);

        const T chiral_shift = 2* M_PI/3;

        T delta_phi = phi_j - phi_i;
        T cos_delta_phi = std::cos(delta_phi);
        T term_phi_diff = 2.0 * m_particle_radius * m_particle_radius * std::pow(1 - cos_delta_phi, 2);

        T cos_j_term = std::cos(phi_j - theta + chiral_shift);
        T cos_i_term = std::cos(phi_i - theta + chiral_shift);
        T term_chiral = 2.0 * m_particle_radius * r * (cos_j_term - cos_i_term);

        // Safety check
        T r_eff_squared = r * r + term_phi_diff + term_chiral;
        // const T min_r_eff_squared = 0.01 * m_sigma * m_sigma;  // Prevent singularities
        // if (r_eff_squared < min_r_eff_squared) {
        //     r_eff_squared = min_r_eff_squared;
        // }
        
        T r_eff = std::sqrt(r_eff_squared);

        // dU/dr_eff
        T dU_dr_eff = 4 * m_epsilon * (
            -12.0 * m_sigma12 * std::pow(r_eff, -13) +
            6.0 * m_sigma6  * std::pow(r_eff, -7)
        );

        // Check for NaN or infinite forces
        if (!std::isfinite(dU_dr_eff)) {
            return {{0, 0}};
        }

        // dr_eff/dr
        T d_r_eff_dr = (r + m_particle_radius * (cos_j_term - cos_i_term)) / r_eff;

        // Corrected dr_eff/dphi_i with proper derivative
        T d_cos_delta_phi = -std::sin(delta_phi);
        T d_cos_i_term = std::sin(phi_i - theta + chiral_shift);

        T d_r_eff_dphi = (
            4.0 * m_particle_radius * m_particle_radius * (1 - cos_delta_phi) * d_cos_delta_phi
            + 2.0 * m_particle_radius * r * d_cos_i_term
        ) / r_eff;

        T force_radial = -dU_dr_eff * d_r_eff_dr;
        T force_angular = -dU_dr_eff * d_r_eff_dphi / r;

        // Cap forces to prevent explosions
        // const T max_force = 1000.0 * m_epsilon / m_sigma;
        // if (std::abs(force_radial) > max_force) {
        //     force_radial = (force_radial > 0 ? 1 : -1) * max_force;
        // }
        // if (std::abs(force_angular) > max_force) {
        //     force_angular = (force_angular > 0 ? 1 : -1) * max_force;
        // }

        // Smooth cutoff near boundary
        // if (r > 0.9 * m_cutoff) {
        //     T cutoff_factor = std::cos(M_PI * (r - 0.9 * m_cutoff) / (0.1 * m_cutoff));
        //     cutoff_factor = 0.5 * (1 + cutoff_factor);
        //     force_radial *= cutoff_factor;
        //     force_angular *= cutoff_factor;
        // }

        return {{force_radial, force_angular}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_particle_radius;
    T m_sigma6, m_sigma12;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class LJChiralPotential {
public:
    LJChiralPotential(T epsilon, T sigma, T cutoff, T particle_radius)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_particle_radius(particle_radius),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)) {}

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T r) const {
        if (r == 0 || r > m_cutoff) return 0;
        
        T phi_i = p1.phi;
        T phi_j = p2.phi;
        T theta = std::atan2(dr[1], dr[0]);

        const T chiral_shift = 2* M_PI/3;

        T delta_phi = phi_j - phi_i;
        T cos_delta_phi = std::cos(delta_phi);
        T term_phi_diff = 2.0 * m_particle_radius * m_particle_radius * std::pow(1 - cos_delta_phi, 2);

        T cos_j_term = std::cos(phi_j - theta + chiral_shift);
        T cos_i_term = std::cos(phi_i - theta + chiral_shift);
        T term_chiral = 2.0 * m_particle_radius * r * (cos_j_term - cos_i_term);

        T r_eff_squared = r * r + term_phi_diff + term_chiral;
        // safety check
        // const T min_r_eff_squared = 0.01 * m_sigma * m_sigma;
        // if (r_eff_squared < min_r_eff_squared) {
        //     r_eff_squared = min_r_eff_squared;
        // }
        
        T r_eff = std::sqrt(r_eff_squared);

        T r6 = std::pow(r_eff, 6);
        T r12 = r6 * r6;
        T LJ = 4 * m_epsilon * (m_sigma12 / r12 - m_sigma6 / r6);

        // Smooth cutoff
        // if (r > 0.9 * m_cutoff) {
        //     T cutoff_factor = std::cos(M_PI * (r - 0.9 * m_cutoff) / (0.1 * m_cutoff));
        //     cutoff_factor = 0.5 * (1 + cutoff_factor);
        //     LJ *= cutoff_factor;
        // }

        return LJ;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_particle_radius;
    T m_sigma6, m_sigma12;
};
