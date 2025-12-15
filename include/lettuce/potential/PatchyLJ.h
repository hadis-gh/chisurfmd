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

template<typename T>
inline T findClosestPatch(const T& phi, const T& gamma, const int& patchNums, const std::vector<T>& patchAngles) {
    T min_theta = std::numeric_limits<T>::max();
    T best_theta = 0;
    for (int i = 0; i < patchNums; ++i) {
        T patchAng = patchAngles[i];
        T theta = wrapAngle(phi + patchAng - gamma);
        T abs_theta = std::abs(theta);
        if (abs_theta < min_theta) {
            min_theta = abs_theta;
            best_theta = theta;
        }
    }
    return best_theta;        
}

// ================  Potential for PatchyLJ ================

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJPotential {
public:
    PatchyLJPotential(T epsilonLJ, T sigmaLJ, T cutoffLJ, T sigmaPatch, int patchNums, 
                      std::string mode)
        : m_epsilonLJ(epsilonLJ), m_sigmaLJ(sigmaLJ), m_cutoffLJ(cutoffLJ), 
          m_sigmaPatch(sigmaPatch), m_patchNums(patchNums), m_mode(mode)
    {
        m_twosigmaPatchSq = 2.0 * sigmaPatch * sigmaPatch;
        m_patchAngles.reserve(m_patchNums);
        for (int a = 0; a < m_patchNums; a++) {
            m_patchAngles.push_back(2 * M_PI * a / m_patchNums);
        }
    }

    T operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoffLJ) return 0;

        T baseLJ = lennardJones(R, m_sigmaLJ, m_epsilonLJ);

        if (R < m_sigmaLJ) {
            return baseLJ;
        }

        const T dx = dr[0];
        const T dy = dr[1];
        const T phi_i = p1.phi;
        const T phi_j = p2.phi;

        T gamma_ij = std::atan2(dy, dx);
        T gamma_ji = wrapAngle(gamma_ij + M_PI);

        T A = 0;
        
        if (m_mode == "closest") {
            T best_theta_i = findClosestPatch(phi_i, gamma_ij, m_patchNums, m_patchAngles);
            T best_theta_j = findClosestPatch(phi_j, gamma_ji, m_patchNums, m_patchAngles);
            
            A = std::exp(-best_theta_i * best_theta_i / m_twosigmaPatchSq) * 
                std::exp(-best_theta_j * best_theta_j / m_twosigmaPatchSq);
        } 
        else { // all patches mode
            T sum_ei = 0;
            for (int l = 0; l < m_patchNums; ++l) {
                T patchAng_i = m_patchAngles[l];
                T theta_i = wrapAngle(phi_i + patchAng_i - gamma_ij);
                sum_ei += std::exp(-theta_i * theta_i / m_twosigmaPatchSq);
            }

            T sum_ej = 0;
            for (int m = 0; m < m_patchNums; ++m) {
                T patch_ang_j = m_patchAngles[m];
                T theta_j = wrapAngle(phi_j + patch_ang_j - gamma_ji);
                sum_ej += std::exp(-theta_j * theta_j / m_twosigmaPatchSq);
            }

            A = sum_ei * sum_ej;
        }

        return baseLJ * A;
    }

private:
    T m_epsilonLJ, m_sigmaLJ, m_cutoffLJ, m_sigmaPatch;
    T m_twosigmaPatchSq;
    int m_patchNums;
    std::string m_mode;
    std::vector<T> m_patchAngles;
};

// ================  Force for PatchyLJ ================

template<typename TParticle, typename T = typename TParticle::value_type>
class PatchyLJForce {
public:
    using value_type = T;

    PatchyLJForce(T epsilonLJ, T sigmaLJ, T cutoffLJ, T sigmaPatch, int patchNums,
                  std::string mode)
        : m_epsilonLJ(epsilonLJ), m_sigmaLJ(sigmaLJ), m_cutoffLJ(cutoffLJ), 
          m_sigmaPatch(sigmaPatch), m_patchNums(patchNums), m_mode(mode)
    {
        m_twosigmaPatchSq = 2.0 * sigmaPatch * sigmaPatch;
        m_patchAngles.reserve(m_patchNums);
        for (int a = 0; a < m_patchNums; a++) {
            m_patchAngles.push_back(2 * M_PI * a / m_patchNums);
        }
    }

    Vec<T, 3> operator()(const TParticle& p1, const TParticle& p2, const Vec<T, 2>& dr, const T R) const {
        if (R <= 0 || R > m_cutoffLJ) return {{0, 0, 0}};

        const T dx = dr[0];
        const T dy = dr[1];
        const T r2 = dx * dx + dy * dy;

        T baseLJ = lennardJones(R, m_sigmaLJ, m_epsilonLJ);
        T dVdR = lennardJonesDerivative(R, m_sigmaLJ, m_epsilonLJ);

        if (R < m_sigmaLJ) {
            T fx = -dVdR * dx / R;
            T fy = -dVdR * dy / R;
            return {{fx, fy, 0}};
        }

        const T phi_i = p1.phi;
        const T phi_j = p2.phi;

        T gamma_ij = std::atan2(dy, dx);
        T gamma_ji = wrapAngle(gamma_ij + M_PI);

        T sum_ei = 0, sum_ej = 0;
        T dsum_ei_dtheta = 0, dsum_ej_dtheta = 0;
        
        if (m_mode == "closest") {
            T best_theta_i = findClosestPatch(phi_i, gamma_ij, m_patchNums, m_patchAngles);
            T best_theta_j = findClosestPatch(phi_j, gamma_ji, m_patchNums, m_patchAngles);
            
            T exp_i = std::exp(-best_theta_i * best_theta_i / m_twosigmaPatchSq);
            sum_ei = exp_i;
            dsum_ei_dtheta = (-2.0 * best_theta_i / m_twosigmaPatchSq) * exp_i;
            
            T exp_j = std::exp(-best_theta_j * best_theta_j / m_twosigmaPatchSq);
            sum_ej = exp_j;
            dsum_ej_dtheta = (-2.0 * best_theta_j / m_twosigmaPatchSq) * exp_j;
        } 
        else { // all patch mode
            for (int l = 0; l < m_patchNums; ++l) {
                T patchAng_i = m_patchAngles[l];
                T theta_i = wrapAngle(phi_i + patchAng_i - gamma_ij);
                T exp_i = std::exp(-theta_i * theta_i / m_twosigmaPatchSq);
                sum_ei += exp_i;
                dsum_ei_dtheta += (-2.0 * theta_i / m_twosigmaPatchSq) * exp_i;
            }

            for (int m = 0; m < m_patchNums; ++m) {
                T patch_ang_j = m_patchAngles[m];
                T theta_j = wrapAngle(phi_j + patch_ang_j - gamma_ji);
                T exp_j = std::exp(-theta_j * theta_j / m_twosigmaPatchSq);
                sum_ej += exp_j;
                dsum_ej_dtheta += (-2.0 * theta_j / m_twosigmaPatchSq) * exp_j;
            }
        }

        // Angular modulation factor and its partial derivatives
        T A = sum_ei * sum_ej;
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

        // Total force components
        T fx = -(dVdR * dRdx * A + baseLJ * dA_dx);
        T fy = -(dVdR * dRdy * A + baseLJ * dA_dy);

        // Torques
        T torque_i = -baseLJ * dA_dtheta_i;  // Torque on particle i
        T torque_j = -baseLJ * dA_dtheta_j;  // Torque on particle j

        return {{fx, fy, torque_j}};
    }

private:
    T m_epsilonLJ, m_sigmaLJ, m_cutoffLJ, m_sigmaPatch;
    T m_twosigmaPatchSq;
    int m_patchNums;
    std::string m_mode;
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
            ("sigmaPatch",   po::value<T>()->default_value(.262),                   "anisotropic strength of potential")
            ("patchNums",    po::value<int>()->default_value(1),                    "Number of patches for geometric LJ")
            ("patchMode",    po::value<std::string>()->default_value("closest"),    "Patch interaction mode: 'all' or 'closest'")
        ;
    }

    static auto force(const po::variables_map &vm) {
        return PatchyLJForce<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["sigmaPatch"].as<T>(),
            vm["patchNums"].as<int>(),
            vm["patchMode"].as<std::string>()
        );
    }

    static auto potential(const po::variables_map &vm) {
        return PatchyLJPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["sigmaPatch"].as<T>(),
            vm["patchNums"].as<int>(),
            vm["patchMode"].as<std::string>()
        );
    }
    
    static auto makePotential(const po::variables_map &vm) {
        return PatchyLJPotential<ParticleOriented<T>>(
            vm["LJepsilon"].as<T>(), 
            vm["LJsigma"].as<T>(), 
            vm["LJcutoff"].as<T>(),
            vm["sigmaPatch"].as<T>(),
            vm["patchNums"].as<int>(),
            vm["patchMode"].as<std::string>()
        );
    }

    using ForceType = PatchyLJForce<ParticleOriented<T>>;
};