import numpy as np
import matplotlib.pyplot as plt

# ____________________Global parameters________________________
EPSILON_LJ = 1.0
SIGMA_LJ = 1.0
CUTOFF_LJ = 3.0

# ____________________Helper Functions________________________

def find_patch_angles(n_patch):
    """Calculate patch angles evenly spaced around a circle"""
    return [2 * np.pi * i / n_patch for i in range(n_patch)]

def wrap_angle(alpha):
    """Wrap angle to range [-π, π]"""
    return (alpha + np.pi) % (2 * np.pi) - np.pi

def closest_patch(phi, gamma, n_patch):
    """Find the closest patch angle"""
    patch_angles = find_patch_angles(n_patch)
    min_theta = np.inf
    best_theta = 0
    
    for patch_ang in patch_angles:
        theta = wrap_angle(phi + patch_ang - gamma)
        abs_theta = np.abs(theta)
        if abs_theta < min_theta:
            min_theta = abs_theta
            best_theta = theta
    
    return best_theta

# ____________________Potential Functions________________________

def lennard_jones_pot(epsilon_LJ, sigma_LJ, r):
    """Calculate Lennard-Jones potential"""
    return 4 * epsilon_LJ * ((sigma_LJ / r) ** 12 - (sigma_LJ / r) ** 6)

def angular_pot(phi_i, phi_j, r_ij, sigma_patch, n_patch, mode='all'):
    """Calculate angular potential between two particles"""
    gamma_ij = 0
    gamma_ji = wrap_angle(gamma_ij + np.pi)
    
    twosigma_patch_sq = 2 * sigma_patch ** 2
    
    v_ang = 0
    if mode == 'closest':
        theta_ij = closest_patch(phi_i, gamma_ij, n_patch)
        theta_ji = closest_patch(phi_j, gamma_ji, n_patch)

        v_ang = (np.exp(-theta_ij ** 2 / twosigma_patch_sq) * 
                np.exp(-theta_ji ** 2 / twosigma_patch_sq))
    elif mode == 'all':
        sum_exp_i = 0
        for patch_i in find_patch_angles(n_patch):
            theta_ij = wrap_angle(phi_i + patch_i - gamma_ij)
            sum_exp_i += np.exp(-theta_ij ** 2 / twosigma_patch_sq)
        
        sum_exp_j = 0
        for patch_j in find_patch_angles(n_patch):
            theta_ji = wrap_angle(phi_j + patch_j - gamma_ji)
            sum_exp_j += np.exp(-theta_ji ** 2 / twosigma_patch_sq)
        
        v_ang = sum_exp_i * sum_exp_j
        
    return v_ang

def total_potential(epsilon_LJ, sigma_LJ, phi_i, phi_j, r_ij, sigma_patch, n_patch, mode='all'):
    """Calculate total potential (Lennard-Jones + Angular)"""
    lj_pot = lennard_jones_pot(epsilon_LJ, sigma_LJ, r_ij)
    ang_pot = angular_pot(phi_i, phi_j, r_ij, sigma_patch, n_patch, mode=mode)
    if r_ij < sigma_LJ:
        return lj_pot
    else:
        return lj_pot * ang_pot

# ____________________Plotting Functions________________________

def plot_angular_comparison(sigma_ang=0.3, n_patch=4, r=1.5):
    """Create comparison plot of angular and total potentials"""
    phi1 = 0.0
    
    phi2_vals = np.linspace(0, 2 * np.pi, 300)
    
    U_vals_all = []
    U_vals_closest = []
    U_vals_total_all = []
    U_vals_total_closest = []
    
    for phi2 in phi2_vals:
        pot_val_all = angular_pot(phi1, phi2, r, sigma_ang, n_patch, mode='all')
        pot_val_closest = angular_pot(phi1, phi2, r, sigma_ang, n_patch, mode='closest')
        pot_val_total_all = total_potential(EPSILON_LJ, SIGMA_LJ, phi1, phi2, r, sigma_ang, n_patch, mode='all')
        pot_val_total_closest = total_potential(EPSILON_LJ, SIGMA_LJ, phi1, phi2, r, sigma_ang, n_patch, mode='closest')
        U_vals_all.append(pot_val_all)
        U_vals_closest.append(pot_val_closest)
        U_vals_total_all.append(pot_val_total_all)
        U_vals_total_closest.append(pot_val_total_closest)
    
    U_vals_all = np.array(U_vals_all)
    U_vals_closest = np.array(U_vals_closest)
    U_vals_total_all = np.array(U_vals_total_all)
    U_vals_total_closest = np.array(U_vals_total_closest)
    
    # Get patch angles for vertical lines
    patch_angles = find_patch_angles(n_patch)
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: Angular potential
    ax_ang = axes[0]
    
    ax_ang.plot(phi2_vals, U_vals_all, color='green', linewidth=2, label='All Patches')
    ax_ang.plot(phi2_vals, U_vals_closest, color='orange', linestyle='--', linewidth=2, label='Closest Patch')
    
    # Add vertical lines at patch positions
    for patch in patch_angles:
        ax_ang.axvline(patch, color='red', linestyle=':', alpha=0.4, linewidth=0.8)
    
    ax_ang.set_xlabel('Particle 2 Orientation (φ₂) [rad]', fontsize=12)
    ax_ang.set_ylabel('Angular Term', fontsize=12)
    ax_ang.set_title(f'Angular Modulation (σ_ang={sigma_ang:.2f}, N={n_patch}, r={r:.2f})', fontsize=14)
    ax_ang.grid(True, alpha=0.3)
    ax_ang.set_xlim(0, 2 * np.pi)
    y_min_ang = min(np.min(U_vals_all), np.min(U_vals_closest))
    y_max_ang = max(np.max(U_vals_all), np.max(U_vals_closest))
    ax_ang.set_ylim(y_min_ang - 0.05, y_max_ang + 0.05)
    ax_ang.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax_ang.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
    ax_ang.legend(loc='best')
    
    # Plot 2: Total potential
    ax_total = axes[1]
    ax_total.plot(phi2_vals, U_vals_total_all, color='royalblue', linewidth=2, label='All Patches')
    ax_total.plot(phi2_vals, U_vals_total_closest, color='coral', linestyle='--', linewidth=2, label='Closest Patch')
    
    # Add horizontal line at y=0
    ax_total.axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=0.5)
    
    ax_total.set_xlabel('Particle 2 Orientation (φ₂) [rad]', fontsize=12)
    ax_total.set_ylabel('Total Potential U(φ₂) [ε]', fontsize=12)
    ax_total.set_title(f'Total Potential vs Orientation (σ_ang={sigma_ang:.2f}, N={n_patch}, r={r:.2f})', fontsize=14)
    ax_total.grid(True, alpha=0.3)
    ax_total.set_xlim(0, 2 * np.pi)
    y_min_total = min(np.min(U_vals_total_all), np.min(U_vals_total_closest))
    y_max_total = max(np.max(U_vals_total_all), np.max(U_vals_total_closest))
    ax_total.set_ylim(y_min_total - 0.1, y_max_total + 0.1)
    ax_total.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax_total.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
    ax_total.legend(loc='best')
    
    plt.tight_layout()
    return fig

def plot_potential_vs_orientation(sigma_ang=0.3, n_patch=4, r=1.5, mode='all'):
    """Create plot of total potential vs orientation"""
    phi1 = 0.0
    
    phi2_vals = np.linspace(0, 2 * np.pi, 300)
    
    U_vals = []
    for phi2 in phi2_vals:
        pot_val = total_potential(EPSILON_LJ, SIGMA_LJ, phi1, phi2, r, sigma_ang, n_patch, mode=mode)
        U_vals.append(pot_val)
    
    U_vals = np.array(U_vals)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(phi2_vals, U_vals, color='royalblue', linewidth=2)
    
    # Add patch position markers
    patch_angles = find_patch_angles(n_patch)
    for patch in patch_angles:
        ax.axvline(patch, color='red', linestyle=':', alpha=0.4, linewidth=0.8)
    
    # Add horizontal line at y=0
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=0.5)
    
    ax.set_xlabel('Particle 2 Orientation (φ₂) [rad]', fontsize=12)
    ax.set_ylabel('Total Potential U(φ₂) [ε]', fontsize=12)
    ax.set_title(f'Potential vs Orientation (σ_ang={sigma_ang:.2f}, N={n_patch}, r={r:.2f}, mode={mode})', 
                 fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 2 * np.pi)
    ax.set_ylim(np.min(U_vals) - 0.1, np.max(U_vals) + 0.1)
    ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
                  ['0', 'π/2', 'π', '3π/2', '2π'])
    plt.tight_layout()
    return fig

def plot_potential_vs_distance(epsilon_LJ=1.0, sigma_LJ=1.0, sigma_ang=0.3, 
                               n_patch=4, phi2=0.0):
    """Create plot of potential vs distance"""
    phi1 = 0.0
    
    r_vals = np.linspace(0.5, 4.0, 300)
    
    U_vals_closest = []
    U_vals_all = []
    U_lj_vals = []
    for r in r_vals:
        pot_val_closest = total_potential(epsilon_LJ, sigma_LJ, phi1, phi2, r, sigma_ang, n_patch, mode='closest')
        pot_val_all = total_potential(epsilon_LJ, sigma_LJ, phi1, phi2, r, sigma_ang, n_patch, mode='all')
        lj_val = lennard_jones_pot(epsilon_LJ, sigma_LJ, r)
        
        U_vals_closest.append(pot_val_closest)
        U_vals_all.append(pot_val_all)
        U_lj_vals.append(lj_val)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(r_vals, U_lj_vals, color='green', linewidth=2, label='Lennard-Jones Potential')
    ax.plot(r_vals, U_vals_all, color='royalblue', linestyle='--', linewidth=2, label='Total Potential (All)')
    ax.plot(r_vals, U_vals_closest, color='coral', linestyle='--', linewidth=2, label='Total Potential (Closest)')
    
    # Add vertical lines for important distances
    ax.axvline(x=sigma_LJ, color='k', linestyle=':', alpha=0.6, linewidth=1.2, label=f'σ_LJ={sigma_LJ}')
    ax.axvline(x=CUTOFF_LJ, color='purple', linestyle=':', alpha=0.6, linewidth=1.2, label=f'Cutoff={CUTOFF_LJ}')
    
    # Add horizontal line at y=0
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=0.5)
    
    ax.set_xlabel('Distance r', fontsize=12)
    ax.set_ylabel('Potential U [ε]', fontsize=12)
    ax.set_title(f'Potential vs Distance (σ_ang={sigma_ang:.2f}, N={n_patch}, φ₂={phi2:.2f})', 
                 fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.5, 4.0)
    y_min = min(np.min(U_vals_all), np.min(U_vals_closest), np.min(U_lj_vals))
    y_max = max(np.max(U_vals_all), np.max(U_vals_closest), np.max(U_lj_vals))
    ax.set_ylim(y_min - 0.5, min(y_max + 0.5, 20))
    ax.legend(loc='best')
    plt.tight_layout()
    return fig

def print_parameters():
    """Print the current simulation parameters"""
    print("Current Patchy Potential Parameters:")
    print(f"  ε_LJ: {EPSILON_LJ}")
    print(f"  σ_LJ: {SIGMA_LJ}")
    print(f"  Cutoff: {CUTOFF_LJ}")
    print("\nVisualization Guide:")
    print("  - Red dotted lines: Patch positions")
    print("  - Black vertical lines: σ_LJ and cutoff distances")
    print("  - All Patches: Sum over all patch combinations")
    print("  - Closest Patch: Only closest patches interact")