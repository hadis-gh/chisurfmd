import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, FloatSlider, IntSlider, Dropdown

# ____________________Global parameters________________________
epsilon_LJ = 1
sigma_LJ = 1
cutoff_LJ = 3.0

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
        
        v_ang = sum_exp_i * sum_exp_j / n_patch**2
        
    return v_ang

def total_potential(epsilon_LJ, sigma_LJ, phi_i, phi_j, r_ij, sigma_patch, n_patch, mode='all'):
    """Calculate total potential (Lennard-Jones + Angular)"""
    lj_pot = lennard_jones_pot(epsilon_LJ, sigma_LJ, r_ij)
    ang_pot = angular_pot(phi_i, phi_j, r_ij, sigma_patch, n_patch, mode=mode)
    if r_ij < sigma_LJ:
        return lj_pot
    else:
        return lj_pot * ang_pot

# ____________________Interactive Plot pure angular ________________________

def plot_angular_potential_scan(sigma_ang_slider=0.3, n_patch_slider=4, r_slider=1.0):
    """Create interactive plot of angular potential and total potential"""
    phi1 = 0.0
    sigma_ang = sigma_ang_slider
    n_patch = n_patch_slider
    r = r_slider
    
    phi2_vals = np.linspace(0, 2 * np.pi, 300)
    
    U_vals_all = []
    U_vals_closest = []
    U_vals_LJ = []
    U_vals_total_all = []
    U_vals_total_closest = []
    
    for phi2 in phi2_vals:
        pot_val_all = angular_pot(phi1, phi2, r, sigma_ang, n_patch, mode='all')
        pot_val_closest = angular_pot(phi1, phi2, r, sigma_ang, n_patch, mode='closest')
        pot_val_total_all = total_potential(epsilon_LJ, sigma_LJ, phi1, phi2, r, sigma_ang, n_patch, mode='all')
        pot_val_total_closest = total_potential(epsilon_LJ, sigma_LJ, phi1, phi2, r, sigma_ang, n_patch, mode='closest')
        pot_val_lj = lennard_jones_pot(epsilon_LJ, sigma_LJ, r)
        
        U_vals_all.append(pot_val_all)
        U_vals_closest.append(pot_val_closest)
        U_vals_total_all.append(pot_val_total_all)
        U_vals_total_closest.append(pot_val_total_closest)
        U_vals_LJ.append(pot_val_lj)
    
    # Convert lists to numpy arrays for easier manipulation
    U_vals_all = np.array(U_vals_all)
    U_vals_closest = np.array(U_vals_closest)
    U_vals_total_all = np.array(U_vals_total_all)
    U_vals_total_closest = np.array(U_vals_total_closest)
    U_vals_LJ = np.array(U_vals_LJ)
    
    patch_angles = find_patch_angles(n_patch)
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    
    # Plot 1: Angular potential
    ax_ang = axes[0]
    
    ax_ang.plot(phi2_vals, U_vals_all, label='All Patches', color='royalblue', linewidth=2)
    ax_ang.plot(phi2_vals, U_vals_closest, label='Closest Patch', color='coral', linestyle='--', linewidth=2)
    
    for patch in patch_angles:
        ax_ang.axvline(patch, color='red', linestyle=':', alpha=0.4, linewidth=0.8)
    
    ax_ang.set_ylabel('Angular Potential', fontsize=12)
    ax_ang.set_title(r'$V_{ang}$'+f'(σ_ang={sigma_ang:.2f}, N={n_patch}, r={r:.2f})', fontsize=14)
    ax_ang.grid(True, alpha=0.3)
    ax_ang.set_xlim(0, 2 * np.pi)
    ax_ang.get_xaxis().set_visible(False)
    
    # Set y-limits
    y_limits = [min(np.min(U_vals_all), np.min(U_vals_closest)) - 0.05, 
                max(np.max(U_vals_all), np.max(U_vals_closest)) + 0.05]
    ax_ang.set_ylim(y_limits)
    
    handles, labels = ax_ang.get_legend_handles_labels()
    unique_labels = set()
    unique_handles = []
    for h, l in zip(handles, labels):
        if l not in unique_labels:
            unique_labels.add(l)
            unique_handles.append(h)
    ax_ang.legend(unique_handles, unique_labels, loc='upper right')
    ax_ang.legend(loc='upper right')

    # Plot 2: Total potential
    ax_total = axes[1]
    ax_total.plot(phi2_vals, U_vals_total_all, label='All Patches', color='royalblue', linewidth=2)
    ax_total.plot(phi2_vals, U_vals_total_closest, label='Closest Patch', color='coral', linestyle='--', linewidth=2)
    ax_total.plot(phi2_vals, U_vals_LJ, label=r'$V_{LJ}$', color='green', linestyle='-', linewidth=3, alpha=.5)
    
    ax_total.axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=0.5)
    
    ax_total.set_xlabel('Particle 2 Rotation angle (φ₂) [rad]', fontsize=12)
    ax_total.set_ylabel('Total Potential', fontsize=12)
    ax_total.set_title(r'$V_{ang}*V_{LJ}$'+f'(σ_ang={sigma_ang:.2f}, N={n_patch}, r={r:.2f})', fontsize=14)
    ax_total.grid(True, alpha=0.3)
    ax_total.set_xlim(0, 2 * np.pi)
    
    # Set y-limits for total potential
    y_limits_total = [min(np.min(U_vals_total_all), np.min(U_vals_total_closest)) - 0.1, 
                      max(np.max(U_vals_total_all), np.max(U_vals_total_closest)) + 0.1]
    ax_total.set_ylim(y_limits_total)
    
    ax_total.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
    ax_total.set_xticklabels(['0', 'π/2', 'π', '3π/2', '2π'])
    ax_total.legend(loc='upper right')
    
    plt.tight_layout()
    plt.show()

# ____________________Interactive Plot Total_phi ________________________

def plot_angular_potential_scan_mode(sigma_ang_slider=0.3, n_patch_slider=4, r_slider=1.0, mode='all'):
    """Create interactive plot of potential"""
    phi1 = 0.0
    sigma_ang = sigma_ang_slider
    n_patch = n_patch_slider
    r = r_slider
    
    phi2_vals = np.linspace(0, 2 * np.pi, 300)
    
    U_vals = []
    for phi2 in phi2_vals:
        pot_val = total_potential(epsilon_LJ, sigma_LJ, phi1, phi2, r, sigma_ang, n_patch, mode=mode)
        U_vals.append(pot_val)
    
    U_vals = np.array(U_vals)
    
    plt.figure(figsize=(10, 6))
    plt.plot(phi2_vals, U_vals, color='royalblue', linewidth=2)
    
    # Add patch position markers
    patch_angles = find_patch_angles(n_patch)
    for patch in patch_angles:
        plt.axvline(patch, color='red', linestyle=':', alpha=0.4, linewidth=0.8)
    
    # Add horizontal line at y=0
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=0.5)
    
    plt.xlabel('Particle 2 Orientation (φ₂) [rad]', fontsize=12)
    plt.ylabel('Total Potential U(φ₂) [ε]', fontsize=12)
    plt.title(f'Potential vs Orientation (σ_ang={sigma_ang:.2f}, N={n_patch}, r={r:.2f}, mode={mode})', 
              fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 2 * np.pi)
    plt.ylim(np.min(U_vals) - 0.1, np.max(U_vals) + 0.1)
    plt.xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi],
               ['0', 'π/2', 'π', '3π/2', '2π'])
    plt.tight_layout()
    plt.show()

# ____________________Interactive Plot Total_r ________________________

def plot_distance_dependence(epsilon_LJ_slider=1.0, sigma_LJ_slider=1.0, sigma_ang_slider=0.3, 
                            n_patch_slider=4, phi2_slider=0.0):
    """Create interactive plot of potential"""
    plt.close('all')
    phi1 = 0.0
    epsilon_LJ = epsilon_LJ_slider
    sigma_LJ = sigma_LJ_slider
    sigma_ang = sigma_ang_slider
    n_patch = n_patch_slider
    phi2 = phi2_slider
    
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
    
    plt.figure(figsize=(10, 6))
    plt.plot(r_vals, U_lj_vals, color='green', linewidth=2, label='Lennard-Jones Potential')
    plt.plot(r_vals, U_vals_all, color='royalblue', linestyle='--', linewidth=2, label='Total Potential (All)')
    plt.plot(r_vals, U_vals_closest, color='coral', linestyle='--', linewidth=2, label='Total Potential (Closest)')
    
    # Add vertical lines for important distances
    plt.axvline(x=sigma_LJ, color='k', linestyle=':', alpha=0.6, linewidth=1.2, label=f'σ_LJ={sigma_LJ}')
    plt.axvline(x=cutoff_LJ, color='purple', linestyle=':', alpha=0.6, linewidth=1.2, label=f'Cutoff={cutoff_LJ}')
    
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.3, linewidth=0.5)
    
    plt.xlabel('Distance r', fontsize=12)
    plt.ylabel('Potential U [ε]', fontsize=12)
    plt.title(f'Potential vs Distance (σ_ang={sigma_ang:.2f}, N={n_patch}, φ₂={phi2:.2f})', 
              fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xlim(0.5, 4.0)
    y_min = min(np.min(U_vals_all), np.min(U_vals_closest), np.min(U_lj_vals))
    y_max = max(np.max(U_vals_all), np.max(U_vals_closest), np.max(U_lj_vals))
    plt.ylim(y_min - 0.5, min(y_max + 0.5, 10))
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()

# ___________________ Comparison angular modulation/ total patchy pot _______________________

def create_interactive_angular_plot():
        interact(plot_angular_potential_scan,
                sigma_ang_slider=FloatSlider(min=0.05, max=1.0, step=0.05, value=0.3,description='σ_ang', continuous_update=True),
                n_patch_slider=IntSlider(min=1, max=12, step=1, value=4, description='# Patches', continuous_update=True),
                r_slider=FloatSlider(min=0.5, max=3.0, step=0.05, value=1.5, description='Distance r', continuous_update=True)
        )

# ____________________________ Mode selection angular factor ________________________________

def create_interactive_angular_plot_mode():
        interact(plot_angular_potential_scan_mode,
                sigma_ang_slider=FloatSlider(min=0.05, max=1.0, step=0.05, value=0.3, description='σ_ang', continuous_update=True),
                n_patch_slider=IntSlider(min=1, max=12, step=1, value=4, description='# Patches', continuous_update=True),
                r_slider=FloatSlider(min=0.5, max=3.0, step=0.05, value=1.5, description='Distance r', continuous_update=True),
                mode=Dropdown(options=['all', 'closest'],value='all', description='Mode')
        )

# _______________________________ Lennard-Jones vs Patchy ___________________________________

def create_interactive_distance_plot():
        interact(plot_distance_dependence,
                epsilon_LJ_slider=FloatSlider(min=0.1, max=5.0, step=0.1, value=1.0, description='ε_LJ', continuous_update=True),    
                sigma_LJ_slider=FloatSlider(min=0.5, max=2.0, step=0.05, value=1.0, description='σ_LJ', continuous_update=True),
                sigma_ang_slider=FloatSlider(min=0.05, max=1.0, step=0.05, value=0.3, description='σ_ang', continuous_update=True),
                n_patch_slider=IntSlider(min=1, max=12, step=1, value=4, description='# Patches', continuous_update=True),
                phi2_slider=FloatSlider(min=0, max=2 * np.pi, step=np.pi/12, value=0.0, description='φ₂', continuous_update=True)
        )