#pragma once

#include <vector>
#include <cmath>
#include "lettuce/Vec.h"
#include "lettuce/ParticleDot.h"
#include "lettuce/ParticleOriented.h"

template<typename TParticle, typename T = typename TParticle::value_type>
class ChiralLJGeometricPotential {
public:
    ChiralLJGeometricPotential(T epsilon, T sigma, T cutoff, T particle_radius)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_particle_radius(particle_radius),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)) {}

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R == 0 || R > m_cutoff) return 0;

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T gamma = std::atan2(dr[1], dr[0]);

        const int patchNum = 3;
        const int N = 3;

        T delta_phi = phi_j - phi_i;
        T cos_delta_phi = std::cos(delta_phi);
        T term_phi_diff = 2.0 * m_particle_radius * m_particle_radius * std::pow(1 - cos_delta_phi, 2);
        
        T total_LJ = 0;

        for (int n = 0; n < patchNum; ++n) {
            T chiral_shift = 2 * M_PI * n /N;

            T cos_j_term = std::cos(phi_j - gamma + chiral_shift);
            T cos_i_term = std::cos(phi_i - gamma + chiral_shift);
            T term_chiral = 2.0 * m_particle_radius * R * (cos_j_term - cos_i_term);

            T r_eff_sq = R * R + term_phi_diff + term_chiral;

            T r3 = std::pow(r_eff_sq, 3);
            T r6 = r3 * r3;

            total_LJ += 4 * m_epsilon * (m_sigma12 / r6 - m_sigma6 / r3);
        }

        return total_LJ;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_particle_radius;
    T m_sigma6, m_sigma12;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class ChiralLJGeometricForce {
public:
    ChiralLJGeometricForce(T epsilon, T sigma, T cutoff, T particle_radius)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), m_particle_radius(particle_radius),
          m_sigma6(std::pow(sigma, 6)), m_sigma12(std::pow(sigma, 12)) {}

    Vec<T, 2> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
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
            T dV_dr_eff = 4 * m_epsilon * (
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

        return {{force_radial, torque}};
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_particle_radius;
    T m_sigma6, m_sigma12;
};
