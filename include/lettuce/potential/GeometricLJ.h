#pragma once

#include <vector>
#include <cmath>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

template<typename T>
T constexpr lennardJones(const T& r, const T& sigma, const T& epsilon) {
    if (r <= 0) return std::numeric_limits<T>::infinity();
    const T sr = sigma / r;
    const T sr6 = sr * sr * sr * sr * sr * sr;
    const T sr12 = sr6 * sr6;
    return 4 * epsilon * (sr12 - sr6);
}

template<typename TParticle, typename T = typename TParticle::value_type>
class GeometricLJPotential {
public:
    GeometricLJPotential(T epsilon, T sigma, T cutoff, T patch_radius, int num_patches = 3)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), 
          m_patch_radius(patch_radius), m_num_patches(num_patches) {

        m_cutoff_sq = cutoff * cutoff;
    }

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return 0;

        T total_potential = lennardJones(R, m_sigma, m_epsilon);

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T rho = m_patch_radius;

        for (int n = 0; n < m_num_patches; ++n) {
            const T theta_i_n = phi_i + 2.0 * M_PI * n / m_num_patches;
            const T theta_j_n = phi_j + 2.0 * M_PI * n / m_num_patches;

            const T delta_cos = rho * (std::cos(theta_j_n) - std::cos(theta_i_n));
            const T delta_sin = rho * (std::sin(theta_j_n) - std::sin(theta_i_n));

            const T r_patch_sq = (dr[0] + delta_cos) * (dr[0] + delta_cos) + 
                                 (dr[1] + delta_sin) * (dr[1] + delta_sin);

            if (r_patch_sq <= m_cutoff_sq) {
                const T r_patch = std::sqrt(r_patch_sq);
                total_potential += lennardJones(r_patch, m_sigma, m_epsilon);
            }
        }

        return total_potential;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_patch_radius;
    T m_cutoff_sq;
    int m_num_patches;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class GeometricLJForce {
public:
    GeometricLJForce(T epsilon, T sigma, T cutoff, T particle_radius)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_particle_radius(particle_radius),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)) {}

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R == 0 || R > m_cutoff) return {{0, 0}};

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T gamma = std::atan2(dr[1], dr[0]);

        const int patchNum = 3;
        const int N = 3;

        T delta_phi = phi_j - phi_i;
        T cos_delta_phi = std::cos(delta_phi);
        T term_phi_diff = 2.0 * m_particle_radius * m_particle_radius * (1 - cos_delta_phi);

        T force_radial = 0;
        T torque = 0;

        for (int n = 0; n < patchNum; ++n) {
            T chiral_shift = n * 2 * M_PI / N;

            T cos_j_term = std::cos(phi_j - gamma + chiral_shift);
            T cos_i_term = std::cos(phi_i - gamma + chiral_shift);
            T term_chiral = 2.0 * m_particle_radius * R * (cos_j_term - cos_i_term);

            T r_eff_sq = R * R + term_phi_diff + term_chiral;

            // dV/dr_eff_sq
            T dV_dr_eff = 4 * m_epsilon / N * (
                -6.0 * m_sigma12 * std::pow(r_eff_sq, -7) +
                 3.0 * m_sigma6  * std::pow(r_eff_sq, -4)
            );

            // dr_eff/dR
            T d_r_eff_dr = 2 * R + 2 * m_particle_radius * (cos_j_term - cos_i_term);

            // dr_eff/dphi_i
            T d_cos_delta_phi = -std::sin(delta_phi);

            T d_r_eff_dphi_i = 4 * m_particle_radius * m_particle_radius * (1 - cos_delta_phi) * d_cos_delta_phi
                             + 2 * m_particle_radius * R * (- std::sin(phi_j - gamma + chiral_shift) + std::sin(phi_i - gamma + chiral_shift));

            force_radial += -dV_dr_eff * d_r_eff_dr;
            torque       += -dV_dr_eff * d_r_eff_dphi_i;
        }

        return {{force_radial*dr[0]/R, force_radial*dr[1]/R, torque}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_particle_radius;
    T m_sigma6, m_sigma12;
};
