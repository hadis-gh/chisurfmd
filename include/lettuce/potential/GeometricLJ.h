#pragma once

#include <vector>
#include <cmath>
#include "lettuce/core/Vec.h"
#include "lettuce/core/ParticleDot.h"
#include "lettuce/core/ParticleOriented.h"

template<typename T>
T constexpr lennardJonesPotential(const T& r, const T& sigma, const T& epsilon) {
    if (r <= 0) return std::numeric_limits<T>::infinity();
    const T sr = sigma / r;
    const T sr6 = sr * sr * sr * sr * sr * sr;
    const T sr12 = sr6 * sr6;
    return 4 * epsilon * (sr12 - sr6);
}

template<typename T>
T constexpr lennardJonesForce(const T& r, const T& sigma, const T& epsilon) {
    if (r <= 0) return std::numeric_limits<T>::infinity();
    const T sr = sigma / r;
    const T sr6 = sr * sr * sr * sr * sr * sr;
    const T sr12 = sr6 * sr6;
    return 24 * epsilon / r * (2 * sr12 - sr6);
}

template<typename TParticle, typename T = typename TParticle::value_type>
class GeometricLJPotential {
public:
    GeometricLJPotential(T epsilon, T sigma, T cutoff, T patch_radius, int num_patches = 3)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), 
          m_patch_radius(patch_radius), m_num_patches(num_patches) {}

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return 0;

        T total_potential = lennardJonesPotential(R, m_sigma, m_epsilon);

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T rho = m_patch_radius;

        for (int n = 0; n < m_num_patches; ++n) {
            const T theta_i_n = phi_i + 2.0 * M_PI * n / m_num_patches;
            const T theta_j_n = phi_j + 2.0 * M_PI * n / m_num_patches;

            const T dx = dr[0] + rho * (std::cos(theta_j_n) - std::cos(theta_i_n));
            const T dy = dr[1] + rho * (std::sin(theta_j_n) - std::sin(theta_i_n));

            const T r_patch_sq = dx * dx + dy * dy;
            const T r_patch = std::sqrt(r_patch_sq);

            if (r_patch <= m_cutoff) {
                total_potential += lennardJonesPotential(r_patch, m_sigma, m_epsilon);
            }
        }

        return total_potential;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_patch_radius;
    int m_num_patches;
};

template<typename TParticle, typename T = typename TParticle::value_type>
class GeometricLJForce {
public:
    GeometricLJForce(T epsilon, T sigma, T cutoff, T patch_radius, int num_patches = 3)
        : m_epsilon(epsilon), m_sigma(sigma), m_cutoff(cutoff), 
          m_patch_radius(patch_radius), m_num_patches(num_patches) {}

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoff) return {{0, 0, 0}};

        Vec<T, 3> f = {{0, 0, 0}};

        // Center-to-center force
        const T f_lj = lennardJonesForce(R, m_sigma, m_epsilon);
        f[0] -= f_lj * dr[0] / R;
        f[1] -= f_lj * dr[1] / R;

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;
        const T rho = m_patch_radius;

        for (int n = 0; n < m_num_patches; ++n) {
            const T theta_i_n = phi_i + 2.0 * M_PI * n / m_num_patches;
            const T theta_j_n = phi_j + 2.0 * M_PI * n / m_num_patches;

            const T dcos = std::cos(theta_j_n) - std::cos(theta_i_n);
            const T dsin = std::sin(theta_j_n) - std::sin(theta_i_n);

            const T dx = dr[0] + rho * dcos;
            const T dy = dr[1] + rho * dsin;
            const T r_patch_sq = dx * dx + dy * dy;
            const T r_patch = std::sqrt(r_patch_sq);

            if (r_patch <= m_cutoff && r_patch > 0) {
                const T f_patch = lennardJonesForce(r_patch, m_sigma, m_epsilon);

                // Force components
                const T fx = -f_patch * dx / r_patch;
                const T fy = -f_patch * dy / r_patch;
                f[0] += fx;
                f[1] += fy;

                // Torques
                const T dtheta_i = rho * (-std::sin(theta_i_n) * fx + std::cos(theta_i_n) * fy);
                const T dtheta_j = rho * ( std::sin(theta_j_n) * fx - std::cos(theta_j_n) * fy);
                f[2] += dtheta_i;
                f[2] += dtheta_j;
            }
        }

        return f;
    }

private:
    T m_epsilon, m_sigma, m_cutoff, m_patch_radius;
    int m_num_patches;
};