import numpy as np
import matplotlib.pyplot as plt
import data_extraction
import system_analysis
import ipywidgets
from ipywidgets import interactive, FloatSlider, IntSlider, Dropdown
import ipywidgets as widgets
from IPython.display import display
#  in the notebook you can try:
# %load_ext autoreload
# %autoreload 2

#------------------------------   POTENTIAL   ------------------------------

SIGMA = 1.0
EPSILON = 1.0

def plot_heatmap_lj_mod_potential(phi_order=1, angular_scale=0.1, angular_order=12):
    """Plot heatmap of potential and force."""
    L = 100
    r = np.linspace(0.9 * SIGMA, 1.3 * SIGMA, L)
    delta_phi = np.linspace(-np.pi, np.pi, L)
    R, Delta_phi = np.meshgrid(r, delta_phi)
    
    potential_values = system_analysis.lj_mod_potential(R, Delta_phi, phi_order, angular_scale, angular_order, SIGMA, EPSILON)
    
    fig, axs = plt.subplots(figsize=(7, 7), dpi=100)
      
    c1 = axs.contourf(R, Delta_phi, potential_values, levels=100, cmap="viridis")
    cbar = fig.colorbar(c1, ax=axs, label='Potential Energy')
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Potential Energy', fontsize=20)
    
    axs.set_title(r'Potential $U(r, \Delta\phi)$', fontsize=24)
    axs.set_xlabel('Distance $r$', fontsize=20)
    axs.set_ylabel('Δφ (rad)', fontsize=20)
    axs.tick_params(axis='both', which='major', labelsize=16)
    
    plt.tight_layout()
    plt.show()

def plot_heatmap_lj_mod_potential_polar(phi_order=1, angular_scale=0.1, angular_order=12):
    """Plot heatmap of potential in polar coordinates"""
    L = 100
    r = np.linspace(0.9 * SIGMA, 1.3 * SIGMA, L)
    delta_phi = np.linspace(0, 2 * np.pi, L)
    R, Delta_phi = np.meshgrid(r, delta_phi)
    
    potential_values = system_analysis.lj_mod_potential(R, Delta_phi, phi_order, angular_scale, angular_order, SIGMA, EPSILON)
    
    fig, axs = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(7, 7))
      
    # Plot the heatmap
    c1 = axs.contourf(Delta_phi, R, potential_values, levels=100, cmap="viridis")
    cbar = fig.colorbar(c1, ax=axs, label='Potential Energy', shrink=0.8)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Potential Energy', fontsize=20)
    
    # Customize the polar plot
    axs.set_title(r'Potential $U(r, \Delta\phi)$', pad=20, fontsize=24)
    axs.set_xlabel('Δφ (rad)', labelpad=20, fontsize=20)
    axs.set_ylabel('', fontsize=20)
    axs.tick_params(axis='both', which='major', labelsize=16)
    
    axs.grid(False)  # Disable the grid
    axs.set_yticklabels([])  # Remove radial tick labels
    
    plt.tight_layout()
    plt.show()

def plot_heatmap_LJmod_pair(potential_type='RRUU', gamma=0, phi_order=1, angular_scale=0.1, angular_order=12):
    L = 100
    R = 0.9 * SIGMA
    phi1 = np.linspace(-np.pi, np.pi, L)
    phi2 = np.linspace(-np.pi, np.pi, L)
    Phi1, Phi2 = np.meshgrid(phi1, phi2)
    
    potential_values = system_analysis.LJmod_pair(R, Phi1, Phi2, phi_order, gamma, potential_type, 
                                      angular_scale, angular_order, SIGMA, EPSILON)
    
    fig, axs = plt.subplots(figsize=(7, 7), dpi=100)
      
    c1 = axs.contourf(Phi1, Phi2, potential_values, levels=100, cmap="viridis")
    # Set shrink to make colorbar shorter
    cbar = fig.colorbar(c1, ax=axs, label='Potential Energy', shrink=0.8)
    cbar.ax.tick_params(labelsize=16)
    cbar.set_label('Potential Energy', fontsize=20)

    axs.set_title(rf'Potential $U(r, \Delta\phi)$, {potential_type}', fontsize=24)
    axs.set_xlabel('φ1 (rad)', fontsize=20)
    axs.set_ylabel('φ2 (rad)', fontsize=20)
    axs.tick_params(axis='both', which='major', labelsize=16)

    axs.set_aspect('equal')
    plt.tight_layout()
    plt.show()

############## interactive plots for modified and chiral Potential ############## 

def create_interactive_plots_LJmod():
    """Create interactive plots with widgets."""
    phi_order_slider = ipywidgets.IntSlider(value=2, min=1, max=10, step=1, description='φ order:')
    angular_scale_slider = ipywidgets.FloatSlider(value=1.6, min=0.1, max=5.0, step=0.1, description='Angular scale:')
    angular_order_slider = ipywidgets.IntSlider(value=12, min=1, max=20, step=1, description='Angular order(r^):')
    
    heatmap_potential = ipywidgets.interactive(
        plot_heatmap_lj_mod_potential,
        phi_order=phi_order_slider,
        angular_scale=angular_scale_slider,
        angular_order=angular_order_slider
    )
    
    heatmap_potential_polar = ipywidgets.interactive(
        plot_heatmap_lj_mod_potential_polar,
        phi_order=phi_order_slider,
        angular_scale=angular_scale_slider,
        angular_order=angular_order_slider
    )
    
    return heatmap_potential, heatmap_potential_polar
    
def create_interactive_plots_LJmod_pair():
    type_slider = ipywidgets.Dropdown(
        options=['RRUU', 'RLUU', 'RRUD', 'RLUD'],
        value='RRUU',
        description='Potential type:'
    )
    gamma_slider = ipywidgets.FloatSlider(
        value=0, 
        min=0, 
        max=2*np.pi, 
        step=np.pi/12, 
        description='γ face Angle:'
    )
    phi_order_slider = ipywidgets.IntSlider(
        value=2, 
        min=1, 
        max=10, 
        step=1, 
        description='φ order:'
    )
    angular_scale_slider = ipywidgets.FloatSlider(
        value=1.6, 
        min=0.1, 
        max=5.0, 
        step=0.1, 
        description='A0:'
    )
    angular_order_slider = ipywidgets.IntSlider(
        value=12, 
        min=1, 
        max=20, 
        step=1, 
        description='m (r^m):'
    )
    
    heatmap_potential = ipywidgets.interactive(
        plot_heatmap_LJmod_pair,
        potential_type=type_slider,
        gamma=gamma_slider,
        phi_order=phi_order_slider,
        angular_scale=angular_scale_slider,
        angular_order=angular_order_slider
    )
        
    return heatmap_potential

#------------------------------   NEW CHIRAL POTENTIAL VISUALIZATION   ------------------------------

def plot_heatmap_lj_chiral_new(
        C=1.0, gamma=0.,
        pot=system_analysis.pair_pot, pot_kwargs={},
        r_range=(0.9 * SIGMA, 1.5 * SIGMA)):
    """Plot heatmap of the new chiral potential U(r, φ2) for fixed φ1=0"""
    L = 100
    r = np.linspace(*r_range, L)
    phi2 = np.linspace(-np.pi, np.pi, L)
    R, Phi2 = np.meshgrid(r, phi2)
    
    phi1_fixed = 0  # fixed orientation for first particle
    potential_values = pot(R, phi1_fixed, Phi2, gamma,
        pot_kwargs={**pot_kwargs, "field":system_analysis.field,
            "field_kwargs":{**pot_kwargs.get('field_kwargs', {}), 'coupling_scale':C}})

    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.contourf(R, Phi2, potential_values, levels=100, cmap='viridis')
    fig.colorbar(c, ax=ax, label='Potential Energy')

    ax.set_title(fr'Chiral Potential $U(r, \varphi_2)$ (C={C}, γ={gamma:.2f})')
    ax.set_xlabel('Distance $r$')
    ax.set_ylabel(r'$\varphi_2$ (rad)')

    plt.tight_layout()
    plt.show()

def plot_heatmap_lj_chiral_new_polar(C=1.0, gamma=0, pot=system_analysis.pair_pot, pot_kwargs={}, r_range=(0.9 * SIGMA, 1.5 * SIGMA)):
    """Plot heatmap of the chiral potential in polar coordinates (φ2 vs r, φ1 fixed)"""
    L = 100
    r = np.linspace(*r_range, L)
    phi2 = np.linspace(0, 2 * np.pi, L)
    R, Phi2 = np.meshgrid(r, phi2)

    phi1_fixed = 0
    potential_values = pot(R, phi1_fixed, Phi2, gamma,
        pot_kwargs={**pot_kwargs, "field":system_analysis.field,
            "field_kwargs":{**pot_kwargs.get('field_kwargs', {}), 'coupling_scale':C}})

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(7, 7))
    c = ax.contourf(Phi2, R, potential_values, levels=100, cmap="viridis")
    fig.colorbar(c, ax=ax, label='Potential Energy')

    ax.set_title(r'Chiral Potential $U(r, \varphi_2)$', pad=20)
    ax.set_yticklabels([])  # Hide radial labels
    ax.grid(False)

    plt.tight_layout()
    plt.show()

def plot_chiral_new(gamma=0, coupling_scale=1.0, coupling_power=6, pot=system_analysis.pair_pot, pot_kwargs={}):
    L = 100
    R = 1.1 * SIGMA
    phi1 = np.linspace(-np.pi, np.pi, L)
    phi2 = np.linspace(-np.pi, np.pi, L)
    Phi1, Phi2 = np.meshgrid(phi1, phi2)
    
    potential_values = pot(R, Phi1, Phi2, gamma,
        pot_kwargs={**pot_kwargs, "field":system_analysis.field,
            "field_kwargs":{**pot_kwargs.get('field_kwargs', {}), 'coupling_scale':coupling_scale, 'coupling_power':coupling_power}})
    
    fig, axs = plt.subplots(figsize=(7, 7), dpi=100)
      
    c1 = axs.contourf(Phi1, Phi2, potential_values, levels=100, cmap="viridis")
    fig.colorbar(c1, ax=axs, label='Potential Energy')
    
    axs.set_title(rf'Potential $U(r, \Delta\phi)$\nType: γ: {gamma:.2f}')
    axs.set_xlabel('φ1 (rad)')
    axs.set_ylabel('φ2 (rad)')
    
    axs.set_aspect('equal')
    plt.tight_layout()
    plt.show()

def plot_angular_psi(coupling_scale=1.0, coupling_power=6, pot=system_analysis.pair_pot_psi, pot_kwargs={}):
    L = 100
    R = 1.1 * SIGMA
    phi1 = np.linspace(-np.pi, np.pi, L)
    phi2 = np.linspace(-np.pi, np.pi, L)
    Phi1, Phi2 = np.meshgrid(phi1, phi2)
    
    if "field_kwargs" not in pot_kwargs:
        pot_kwargs['field_kwargs'] = {}
    pot_kwargs['field_kwargs']['coupling_scale'] = coupling_scale
    pot_kwargs['field_kwargs']['coupling_power'] = coupling_power

    potential_values = pot(R, Phi1, Phi2,
        **pot_kwargs)
    
    fig, axs = plt.subplots(figsize=(7, 7), dpi=100)
      
    c1 = axs.contourf(Phi1, Phi2, potential_values, levels=100, cmap="viridis")
    fig.colorbar(c1, ax=axs, label='Potential Energy')
    
    axs.set_title(rf'Potential $U(r, \Delta\phi)$')
    axs.set_xlabel(rf'$\psi_1$ (rad)')
    axs.set_ylabel(rf'$\psi_2$ (rad)')

    axs.set_aspect('equal')
    plt.tight_layout()
    plt.show()


# def plot_chiral_new_pair(potential_type='RRUU', gamma=0, coupling_scale=1.0, coupling_power=6):
#     L = 100
#     R = 1.1 * SIGMA
#     phi1 = np.linspace(-np.pi, np.pi, L)
#     phi2 = np.linspace(-np.pi, np.pi, L)
#     Phi1, Phi2 = np.meshgrid(phi1, phi2)
    
#     potential_values = system_analysis.chiral_new_pair(R, Phi1, Phi2, gamma, potential_type, coupling_scale, coupling_power, SIGMA, EPSILON)
    
#     fig, axs = plt.subplots(figsize=(7, 7), dpi=100)
      
#     c1 = axs.contourf(Phi1, Phi2, potential_values, levels=100, cmap="viridis")
#     fig.colorbar(c1, ax=axs, label='Potential Energy')
    
#     axs.set_title(rf'Potential $U(r, \Delta\phi)$\nType: {potential_type}, γ: {gamma:.2f}')
#     axs.set_xlabel('φ1 (rad)')
#     axs.set_ylabel('φ2 (rad)')
    
#     axs.set_aspect('equal')
#     plt.tight_layout()
#     plt.show()

############## interactive plots for new chiral Potential ############## 

def create_interactive_chiral_new(pot=system_analysis.pair_pot, pot_kwargs={}):
    """Interactive widgets for the new chiral potential"""
    C_slider = ipywidgets.FloatSlider(value=10.0, min=0.1, max=50.0, step=0.1, description='C:')
    gamma_slider = ipywidgets.FloatSlider(value=0, min=0, max=2*np.pi, step=np.pi/12, description='γ:')

    heatmap = ipywidgets.interactive(
        lambda C, gamma: plot_heatmap_lj_chiral_new(C, gamma, pot=pot, pot_kwargs=pot_kwargs),
        C=C_slider,
        gamma=gamma_slider,
    )

    polar_plot = ipywidgets.interactive(
        lambda C, gamma: plot_heatmap_lj_chiral_new_polar(C, gamma, pot=pot, pot_kwargs=pot_kwargs),
        C=C_slider,
        gamma=gamma_slider
    )

    return heatmap, polar_plot

def create_interactive_plots_chiral_new_pair(pot=system_analysis.pair_pot, pot_kwargs={}):
    type_slider = ipywidgets.Dropdown(
        options=['RRUU', 'RLUU', 'RRUD', 'RLUD'],
        value='RRUU',
        description='Potential type:'
    )
    gamma_slider = ipywidgets.FloatSlider(
        value=0, 
        min=0, 
        max=2*np.pi, 
        step=np.pi/12, 
        description='γ face Angle:'
    )
    coupling_scale_slider = ipywidgets.FloatSlider(
        value=1.6, 
        min=0.1, 
        max=5.0, 
        step=0.1, 
        description='coupling scale:'
    )
    coupling_order_slider = ipywidgets.IntSlider(
        value=12, 
        min=1, 
        max=20, 
        step=1, 
        description='m (r^m):'
    )

    heatmap_potential = ipywidgets.interactive(
        lambda gamma, coupling_scale, coupling_order: plot_chiral_new(gamma, coupling_scale=coupling_scale, coupling_power=coupling_order, pot=pot, pot_kwargs=pot_kwargs),
        gamma=gamma_slider,
        coupling_scale=coupling_scale_slider,
        coupling_order=coupling_order_slider
    )
        
    return heatmap_potential

def create_interactive_plots_chiral_new_pair_psi(pot=system_analysis.pair_pot_psi, pot_kwargs={}):
    coupling_scale_slider = ipywidgets.FloatSlider(
        value=1.6, 
        min=0.1, 
        max=5.0, 
        step=0.1, 
        description='coupling scale:'
    )
    coupling_order_slider = ipywidgets.IntSlider(
        value=12, 
        min=1, 
        max=20, 
        step=1, 
        description='m (r^m):'
    )

    heatmap_potential = ipywidgets.interactive(
        lambda coupling_scale, coupling_order: plot_angular_psi(coupling_scale=coupling_scale, coupling_power=coupling_order, pot=pot, pot_kwargs=pot_kwargs),
        coupling_scale=coupling_scale_slider,
        coupling_order=coupling_order_slider
    )
        
    return heatmap_potential


#------------------------------   SINGLE MD   ------------------------------

############## plot system properties, Energy, Temperature, Neighbors, Order, COM vel ##############

def plot_temperature(m_temperatures, temperature_label=False):
    temperature_data = m_temperatures['data']
    temp_labels = m_temperatures['temperature label']
    fig, ax = plt.subplots(figsize=(7, 4))
    
    ax.grid(True, linestyle='--', alpha=0.7, color='gray')

    colors = plt.cm.plasma(np.linspace(0, 1, 4))
    axis_labels = ['x', 'y', r'$\omega$']

    for i, axis_label in enumerate(axis_labels):
        ax.plot(temperature_data[:, :, i], label=f"T{axis_label}", color=colors[i])
    
    if temperature_label: 
        total_points = temperature_data.shape[0]
        tick_indices = np.linspace(0, total_points - 1, 10, dtype=int)
        label_indices = np.linspace(0, len(temp_labels) - 1, 10, dtype=int)
        ax.set_xlabel("Temperature")
        ax.set_xticks(tick_indices)
        ax.set_xticklabels([f"T={temp_labels[i][0]:.1f}" for i in label_indices], rotation=45)
    else:
        ax.set_xlabel("Steps")

                
    ax.set_ylabel("Measured Temperature")
    ax.set_title("Temperature Evolution")
    
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=True, fancybox=True, shadow=False)
    
    plt.tight_layout()
    plt.show()

def plot_neighbors(neighbors_array, m_temperatures, temperature_show=False, temperature_label=False):
    neighbors_data = neighbors_array['data']
    temp_labels = neighbors_array['temperature label']
    fig, ax = plt.subplots(figsize=(7, 4))

    ax.grid(True, linestyle='--', alpha=0.7, color='gray')

    colors = plt.cm.viridis(np.linspace(0, 1, 4))

    for i in range(neighbors_data.shape[2]):
        ax.plot(neighbors_data[:, :, i], label=f"Shell {i+1}", color=colors[i], linewidth=2)


    if temperature_label: 
        total_points = neighbors_data.shape[0]
        tick_indices = np.linspace(0, total_points - 1, 10, dtype=int)
        label_indices = np.linspace(0, len(temp_labels) - 1, 10, dtype=int)
        ax.set_xlabel("Temperature")
        ax.set_xticks(tick_indices)
        ax.set_xticklabels([f"T={temp_labels[i][0]:.1f}" for i in label_indices], rotation=45)
    else:
        ax.set_xlabel("Steps")

    ax.set_ylabel("Number of Neighbors")
    ax.set_title("Neighbor Analysis")
    
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=True, fancybox=True, shadow=False)

    plt.tight_layout()
    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

def plot_energies(kinetic_energy, potential_energy, m_temperatures, show_potential=True, temperature_show=False, temperature_label=False, plot_window=None):
    if plot_window is None:     
        plot_window = kinetic_energy['data'].shape[0]
    kinetic_energy_data = kinetic_energy['data'][:plot_window]
    temp_labels = kinetic_energy['temperature label'][:plot_window]
    potential_energy_data = potential_energy['data'][:plot_window]

    # kinetic_energy_sum = np.sum(kinetic_energy_data, axis=2)[:, 0]
    # if y2.shape[0] > kinetic_energy_sum.shape[0]:
    #     y2 = y2[:kinetic_energy_sum.shape[0]]  # Truncate y2
    # else:
    #     kinetic_energy_sum = kinetic_energy_sum[:y2.shape[0]]  # Truncate the sum

    # ax.plot(y2 + kinetic_energy_sum, label='Total', color=colors[4])

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.grid(True, linestyle='--', alpha=0.7, color='gray')

    colors = plt.cm.plasma(np.linspace(0, 1, kinetic_energy_data.shape[2] + 3))
    axis_labels = ['x', 'y', r'$\omega$']

    for i, axis_label in enumerate(axis_labels):
        ax.plot(kinetic_energy_data[:, :, i], label=f"K{axis_label}", color=colors[i])

    if show_potential:
        y2 = potential_energy_data.flatten()
        ax.plot(y2, label='U', color=colors[3])
        ax.plot(y2 + np.sum(kinetic_energy_data, axis=2)[:, 0], label='Total', color=colors[4])

    if temperature_label: 
        total_points = kinetic_energy_data.shape[0]
        tick_indices = np.linspace(0, total_points - 1, 10, dtype=int)
        label_indices = np.linspace(0, len(temp_labels) - 1, 10, dtype=int)
        ax.set_xlabel("Temperature")
        ax.set_xticks(tick_indices)
        ax.set_xticklabels([f"T={temp_labels[i][0]:.1f}" for i in label_indices], rotation=45)
    else:
        ax.set_xlabel("Steps")

    ax.set_ylabel("Energy")
    ax.set_title("Energy Evolution")
    
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=True, fancybox=True, shadow=False)

    plt.tight_layout()
    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

def plot_com_velocity(com_velocity, m_temperatures, temperature_show=False, temperature_label=False):
    com_vel_data = com_velocity['data']
    temp_labels = com_velocity['temperature label']
    fig, ax = plt.subplots(figsize=(7, 4))

    ax.grid(True, linestyle='--', alpha=0.7, color='gray')

    colors = plt.cm.plasma(np.linspace(0, 1, 4))
    axis_labels = ['x', 'y']

    for i, axis_label in enumerate(axis_labels):
        ax.plot(com_vel_data[:, :, i], label=f"K{axis_label}", color=colors[i])

    if temperature_label: 
        total_points = com_vel_data.shape[0]
        tick_indices = np.linspace(0, total_points - 1, 10, dtype=int)
        label_indices = np.linspace(0, len(temp_labels) - 1, 10, dtype=int)
        ax.set_xlabel("Temperature")
        ax.set_xticks(tick_indices)
        ax.set_xticklabels([f"T={temp_labels[i][0]:.1f}" for i in label_indices], rotation=45)
    else:
        ax.set_xlabel("Steps")

    ax.set_ylabel("Velocity")
    ax.set_title("Center of Mass Velocity Evolution")
    
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left', frameon=True, fancybox=True, shadow=False)

    plt.tight_layout()
    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

def plot_com_ang_velocity(com_ang_velocity, m_temperatures, temperature_show=False, temperature_label=False):
    data = com_ang_velocity['data']
    temp_labels = com_ang_velocity['temperature label']

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.grid(True, linestyle='--', alpha=0.7, color='gray')

    y2 = data.flatten()
    ax.plot(y2, color='#b47c39')

    if temperature_label: 
        total_points = data.shape[0]
        tick_indices = np.linspace(0, total_points - 1, 10, dtype=int)
        label_indices = np.linspace(0, len(temp_labels) - 1, 10, dtype=int)
        ax.set_xlabel("Temperature")
        ax.set_xticks(tick_indices)
        ax.set_xticklabels([f"T={temp_labels[i][0]:.1f}" for i in label_indices], rotation=45)
    else:
        ax.set_xlabel("Steps")

    ax.set_ylabel("Angular Velocity")
    ax.set_title("Center of Mass Angular Velocity Evolution")
    
    plt.tight_layout()
    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)
        
def plot_order_parameter(order, m_temperatures, temperature_show=False, temperature_label=False):
    order_data = order['data']
    temp_labels = order['temperature label']
    fig, ax = plt.subplots(figsize=(6, 3.8))  
    
    ax.grid(True, linestyle='--', alpha=0.7, color='gray')
    ax.plot(order_data[:, :], color='#5499C7')

    if temperature_label: 
        total_points = order_data.shape[0]
        tick_indices = np.linspace(0, total_points - 1, 10, dtype=int)
        label_indices = np.linspace(0, len(temp_labels) - 1, 10, dtype=int)
        ax.set_xlabel("Temperature")
        ax.set_xticks(tick_indices)
        ax.set_xticklabels([f"T={temp_labels[i][0]:.1f}" for i in label_indices], rotation=45)
    else:
        ax.set_xlabel("Steps")

    ax.set_ylabel("Orientation Order")
    ax.set_title(f"Orientation Order Parameter")

    plt.tight_layout()
    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

############## snapshot of system over steps/ single file --one MD simulation or one aggregation ##############

def plot_configuration(positions, handedness=None, step_target=-1, radius=0.4, color_phi=False, diameter=20, 
                       x_limit=(0, 50),y_limit=(0, 50), aggregation=False, 
                       save_fig=False, save_name='output.pdf', 
                       set_title=None):

    x0, x1 = x_limit
    y0, y1 = y_limit
        
    positions_data = positions['data']
    handedness_data = handedness['data'] if handedness is not None else None

    availble_steps = positions_data.shape[0]
    print('available steps: ', availble_steps, '\tchoosed step: ', step_target)

    x = positions_data[step_target, :, 0]
    y = positions_data[step_target, :, 1]
    h = handedness_data[step_target, :] if handedness_data is not None else None

    if positions_data.shape[2] == 3:
        phi = positions_data[step_target, :, 2]

    T = positions["temperature label"][0][0]
    
    fig, ax = plt.subplots(figsize=(6, 4), dpi=150)

    if positions_data.shape[2] == 3:
        if color_phi:
            scatter = ax.scatter(
                x, y, 
                s=diameter,
                c=phi,
                cmap='twilight',
                alpha=0.8
            )
            cbar = fig.colorbar(scatter, ax=ax, label='Phi (radians)', orientation='vertical', shrink=0.8, pad=0.02)
            cbar.set_alpha(1)
        else:
            if aggregation:
                colors = np.arange(len(x))
                colors[:20] = 0
                
                scatter = ax.scatter(
                    x, y, 
                    s=diameter,
                    edgecolors='black',
                    facecolors='none',
                    alpha=0.8,
                    cmap='plasma',
                    c=colors
                )
                x_end = x + radius * np.cos(phi)
                y_end = y + radius * np.sin(phi)
                cbar = fig.colorbar(scatter, ax=ax, label='deposited particles order', orientation='vertical', shrink=0.8, pad=0.02)
                
            else:
                colors_h = ['lightcoral' if val == 1 else 'lightsteelblue' for val in h]
                scatter = ax.scatter(
                    x, y, 
                    s=diameter,
                    edgecolors='black',
                    facecolors= colors_h,
                    alpha=0.8
                )
                x_end = x + radius * np.cos(phi)
                y_end = y + radius * np.sin(phi)
                
            for i in range(positions_data.shape[1]):
                ax.plot([x[i], x_end[i]], [y[i], y_end[i]], color='black', linewidth=0.8, alpha=0.8)

    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    
    if set_title:
        ax.set_title(set_title)
    else:
        if step_target == -1:
            ax.set_title(f"Particle's configuration - Last step")
        elif step_target == 0:
            ax.set_title(f"Particle's configuration - First step")
        else:
            ax.set_title(f"Particle's configuration - {step_target/availble_steps * 100:.0f}%")

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)

    ax.set_aspect('equal', adjustable='box')
    ax.grid(linestyle='--', alpha=0.5)
    ax.set_axisbelow(True)

    plt.tight_layout()
    if save_fig:
        plt.savefig(save_name)
    plt.show()
    
############## animation of system configuration evolution over time steps for single simple file ##############

import cv2

def draw_grid(frame, frame_size, area, grid_spacing=10, color=(200, 200, 200)):
    scale = frame_size / area
    spacing_pixels = int(grid_spacing * scale)

    for x in range(0, frame_size, spacing_pixels):
        cv2.line(frame, (x, 0), (x, frame_size), color, 1)

    for y in range(0, frame_size, spacing_pixels):
        cv2.line(frame, (0, y), (frame_size, y), color, 1)

def map_to_frame(x, y, frame_size, area):
    scale = frame_size / area
    return int(x * scale), int(y * scale)

def animate_position_simple(positions, output_name, line_length=0.4, area=50, color_p=False, frame_skip=1, frame_size=800, fps=20):
    positions_data = positions['data']
    num_steps = positions_data.shape[0]
    num_particles = positions_data.shape[1]
    has_phi = positions_data.shape[2] == 3

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(output_name, fourcc, fps, (frame_size, frame_size))

    particle_radius = max(5, frame_size // (2 * area))

    for step in range(0, num_steps, frame_skip):
        frame = np.ones((frame_size, frame_size, 3), dtype=np.uint8) * 255  # White background
        
        draw_grid(frame, frame_size, area, grid_spacing=10, color=(200, 200, 200))

        for i in range(num_particles):
            x, y = positions_data[step, i, :2]
            
            # Skip NaN values (particle not yet deposited)
            if np.isnan(x) or np.isnan(y):
                continue

            cx, cy = map_to_frame(x, y, frame_size, area)

            if has_phi:
                phi = positions_data[step, i, 2]
            
            # Draw circle with black outline and white fill
            cv2.circle(frame, (cx, cy), particle_radius, (0, 0, 0), 2, lineType=cv2.LINE_AA)
            cv2.circle(frame, (cx, cy), particle_radius - 1, (255, 255, 255), -1, lineType=cv2.LINE_AA)

            if has_phi and not color_p and not np.isnan(phi):
                # Draw direction line from center to end of radius
                x_end = x + line_length * np.cos(phi)
                y_end = y + line_length * np.sin(phi)

                # Skip NaN values in phi-based calculations
                if np.isnan(x_end) or np.isnan(y_end):
                    continue

                cx_end, cy_end = map_to_frame(x_end, y_end, frame_size, area)
                cv2.line(frame, (cx, cy), (cx_end, cy_end), (0, 0, 0), 2, lineType=cv2.LINE_AA)

        # Add text for step count
        text = f"Step {step}/{num_steps}"
        cv2.putText(frame, text, (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 0, 0), 2, lineType=cv2.LINE_AA)

        out.write(frame)

        if step % (50 * frame_skip) == 0:
            print(f"Processing frame{step}/{num_steps}      ", end='\r')

    out.release()
    print("Processing complete.            ", end='\r')
    print("Video saved as", output_name)

def animate_position_handedness(positions, handedness, output_name, line_length=0.38, area=20, color_p=False, frame_skip=1000, frame_size=800, fps=20):
    positions_data = positions['data']
    num_steps = positions_data.shape[0]
    num_particles = positions_data.shape[1]
    has_phi = positions_data.shape[2] == 3

    handedness_data = handedness['data']

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(output_name, fourcc, fps, (frame_size, frame_size))

    particle_radius = max(5, frame_size // int(2.5* area))
    for step in range(0, num_steps, frame_skip):
        frame = np.ones((frame_size, frame_size, 3), dtype=np.uint8) * 255  # White background
        
        draw_grid(frame, frame_size, area, grid_spacing=10, color=(200, 200, 200))

        for i in range(num_particles):
            x, y = positions_data[step, i, :2]
            h = handedness_data[step, i]
            
            color = [(255,204,204) if h == 1 else (205,230,255)]

            # Skip NaN values (particle not yet deposited)
            if np.isnan(x) or np.isnan(y):
                continue

            cx, cy = map_to_frame(x, y, frame_size, area)

            if has_phi:
                phi = positions_data[step, i, 2]
            
            # Draw circle with black outline and white fill
            cv2.circle(frame, (cx, cy), particle_radius, (0, 0, 0), 2, lineType=cv2.LINE_AA)
            cv2.circle(frame, (cx, cy), particle_radius - 1, color[0], -1, lineType=cv2.LINE_AA)

            if has_phi and not color_p and not np.isnan(phi):
                # Draw direction line from center to end of radius
                x_end = x + line_length * np.cos(phi)
                y_end = y + line_length * np.sin(phi)

                # Skip NaN values in phi-based calculations
                if np.isnan(x_end) or np.isnan(y_end):
                    continue

                cx_end, cy_end = map_to_frame(x_end, y_end, frame_size, area)
                cv2.line(frame, (cx, cy), (cx_end, cy_end), (0, 0, 0), 2, lineType=cv2.LINE_AA)

        # Add text for step count
        text = f"Step {step}/{num_steps}"
        cv2.putText(frame, text, (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.8, (0, 0, 0), 2, lineType=cv2.LINE_AA)

        out.write(frame)

        if step % (50 * frame_skip) == 0:
            print(f"Processing frame{step}/{num_steps}      ", end='\r')

    out.release()
    print("Processing complete.            ", end='\r')
    print("Video saved as", output_name)

############## Histogram of φ and Δφ at specific step ##############

def plot_trajectory_steps(positions, step_target, step_window=10000, area=50):
    positions_data = positions['data']
    print('available steps: ', positions_data.shape[0])
    
    step_max = step_target + step_window

    fig, ax = plt.subplots(1, 2, figsize=(8, 4), dpi=150)

    num_particles = positions_data.shape[1]

    for i in range(num_particles):
        ax[0].scatter(
            positions_data[step_target:step_max, i, 0], 
            positions_data[step_target:step_max, i, 1], 
            c=np.linspace(0, 1, step_max - step_target),
            cmap="plasma",
            s=0.02, alpha=0.6
        )

    ax[0].set_title(f"Trajectories, for step: {step_target} to {step_max}")
    ax[0].set_xlabel("X Position")
    ax[0].set_ylabel("Y Position")

    # Plot initial vs final positions
    ax[1].scatter(
        positions_data[step_target, :, 0], positions_data[step_target, :, 1], 
        marker='o', s=10, color='#00aabb', label='Initial Position', alpha=0.8
    )
    ax[1].scatter(
        positions_data[step_max - 1, :, 0], positions_data[step_max - 1, :, 1], 
        marker='x', s=10, color='#ff7777', label='Final Position', alpha=0.8
    )
    print("available steps: ", positions_data.shape[0])
    ax[1].set_title("Initial vs. Final Positions")
    ax[1].set_xlabel("X Position")
    ax[1].set_ylabel("Y Position")
    ax[1].legend()

    for ax_i in ax:
        ax_i.set_xlim(0,area)
        ax_i.set_ylim(0,area)
        ax_i.set_aspect('equal')
        ax_i.grid(linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.show()

def plot_hist_phi_steps(positions, step_target, step_window=1000, bins_num=100, dpi=120):
    positions_data = positions['data']
    
    step_max = step_target + step_window
    print(f"available steps: {positions_data.shape[0]}, chosen step: {step_target} - {step_max}")
    
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize=(6, 4), dpi=dpi)

    ax.hist(positions_data[step_target:step_max,:,2].flatten(), bins=bins_num, color='#954965', alpha=0.7)

    ax.set_title(r'Histogram of particles $\phi$')
    ax.set_xlabel(r'Rotational Angle $\phi$')
    # ax.set_ylabel('Frequency')

    ax.set_axisbelow(True)
    ax.grid(color='gray', linestyle='dashed', alpha=0.5)

    ax.set_theta_zero_location("E")

    plt.tight_layout()
    plt.show()

def plot_deltaphi_hist_steps(positions, step_target, min_dis=20, step_window=1000, bins_num= 100, dpi=120):
    
    step_max = step_target + step_window
    print(f"available steps: {positions['data'].shape[0]}, chosen step: {step_target} - {step_max}")
    
    total_delta_phi=[]
    
    for a in range(step_window):
        data = positions['data'][step_target+a]
        for i in range(data.shape[0]):
            for j in range(i+1, data.shape[0]):
                r = np.sqrt((data[i, 0] - data[j, 0])**2 + (data[i, 1] - data[j, 1])**2)
                if r < min_dis:
                    delta_phi = data[i, 2] - data[j, 2]
                    delta_phi = (delta_phi + np.pi) % (2 * np.pi) - np.pi
                    total_delta_phi.append(delta_phi)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize=(6, 4), dpi=dpi)
    
    bins = np.linspace(-np.pi, np.pi, bins_num + 1)
    hist, bin_edges = np.histogram(total_delta_phi, bins=bins)

    ax.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), color='#55b3d1', align='edge', alpha=0.7)

    ax.set_title(r'Histogram of $\Delta \phi$')
    ax.set_xlabel(r'Rotational Angle $\phi$')
    
    ax.set_axisbelow(True)
    ax.grid(color='gray', linestyle='dashed', alpha=0.5)
    
    ax.set_theta_zero_location("E")
    
    plt.tight_layout()
    plt.show()

#------------------------------   TEMPERATURE LOOP   ------------------------------

############## snapshot of configuration at specific Temperature ##############

def target_temperature_index (variable_name, target_temperature):
    data = variable_name['data']
    temperature_labels = variable_name['temperature label'].flatten()
    print(f"available temperatures are from {min(temperature_labels)} to {max(temperature_labels)}.")
          
    try:
        index = np.where(temperature_labels == target_temperature)[0][0]
    except IndexError:
        print(f"🞩🞩🞩 Target temperature {target_temperature} not found in temperature_labels. Please choose a valid temperature! 🞩🞩🞩")
        return 0

    total_snapshots = data.shape[0]
    return int(total_snapshots / len(temperature_labels) * index)

def plot_trajectory_temperature(positions, target_temperature, step_window=10000):
    positions_data = positions['data']

    step_min = target_temperature_index(positions, target_temperature)
    step_max = step_min + step_window

    fig, ax = plt.subplots(1, 2, figsize=(8, 4), dpi=200)

    num_particles = positions_data.shape[1]

    for i in range(num_particles):
        ax[0].scatter(
            positions_data[step_min:step_max, i, 0], 
            positions_data[step_min:step_max, i, 1], 
            c=np.linspace(0, 1, step_max - step_min),
            cmap="plasma",
            s=0.02, alpha=0.6
        )

    ax[0].set_title(f"Trajectories at T={target_temperature} for {step_window} steps")
    ax[0].set_xlabel("X Position")
    ax[0].set_ylabel("Y Position")

    # Plot initial vs final positions
    ax[1].scatter(
        positions_data[step_min, :, 0], positions_data[step_min, :, 1], 
        marker='o', s=20, color='#00aabb', label='Initial Position'
    )
    ax[1].scatter(
        positions_data[step_max - 1, :, 0], positions_data[step_max - 1, :, 1], 
        marker='x', s=20, color='#ff7777', label='Final Position'
    )
    print("available steps: ", positions_data.shape[0])
    ax[1].set_title("Initial vs. Final Positions")
    ax[1].set_xlabel("X Position")
    ax[1].set_ylabel("Y Position")
    ax[1].legend()

    for ax_i in ax:
        ax_i.set_xlim(0,20)
        ax_i.set_ylim(0,20)
        ax_i.set_aspect('equal', adjustable='datalim')
        ax_i.grid(linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.show()

import numpy as np
import matplotlib.pyplot as plt

def plot_snapshot_temperature(positions, handedness, target_temperature, patchAngs=None,
                               line_length=0.45, area=20, radius=True, color_palette='hsv', dpi=120):
    positions_data = positions['data']
    handedness_data = handedness['data']
    shot = target_temperature_index(positions, target_temperature)

    fig, ax = plt.subplots(figsize=(6, 4), dpi=dpi)
    
    x = positions_data[shot, :, 0]
    y = positions_data[shot, :, 1]
    phi = positions_data[shot, :, 2]  # orientation of each particle
    h = handedness_data[shot, :]

    # Define colors based on handedness
    colors = ['lightcoral' if val == 1 else 'lightsteelblue' for val in h]
    
    # Plot main particles
    ax.scatter(
        x, y, 
        s=110,
        edgecolors='black',
        facecolor=colors,
        alpha=0.8
    )

    if radius:
        # Plot orientation lines
        x_end = x + line_length * np.cos(phi)
        y_end = y + line_length * np.sin(phi)
        for i in range(len(x)):
            ax.plot([x[i], x_end[i]], [y[i], y_end[i]], color='black', linewidth=0.8, alpha=0.8)
    
    elif patchAngs is not None:
        # Plot patches as tiny points
        for i in range(len(x)):
            if np.isnan(phi[i]):
                continue
            for ang in patchAngs[i]:  # support multiple patch angles per particle
                patch_angle = phi[i] + ang  # global angle
                patch_x = x[i] + line_length * np.cos(patch_angle)
                patch_y = y[i] + line_length * np.sin(patch_angle)
                ax.scatter(patch_x, patch_y, s=5, color='blue', alpha=0.9, zorder=5)
    
    else:
        # Fallback: color by orientation
        scatter = ax.scatter(
            x, y,
            c=phi,
            s=60,
            alpha=0.8,
            vmin=-np.pi,
            vmax=np.pi,
            cmap=color_palette
        )
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('φ in radian')

    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    ax.set_title(f'Configuration (T={target_temperature})')
    
    ax.set_xlim(0, area)
    ax.set_ylim(0, area)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(linestyle='--', alpha=0.5)
    ax.set_axisbelow(True)
    plt.tight_layout()
    plt.show()

############## histogram of φ and Δφ at specific Temperature ##############

def plot_hist_phi_temperature(positions, target_t, step_window, bins_num=100, dpi=120):
    step_target = target_temperature_index(positions, target_t)
    plot_hist_phi_steps(positions, step_target, step_window, bins_num, dpi)

def plot_delta_phi_hist_temperature(positions, target_t, min_dis, step_window, bins_num=100, dpi=120):
    step_target = target_temperature_index(positions, target_t)
    plot_deltaphi_hist_steps(positions, step_target, min_dis, step_window, bins_num, dpi)

#------------------------------   AGGREGATION   ------------------------------

############## snapshot of system over steps/ aggregation with specific Temperature & Deposition ##############

def plot_configuration_aggregation(output_dir, temperature, deposition_rate, step_target=-1, 
                                   radius=0.4, diameter=20, x_limit=(0, 50), y_limit=(0, 50), 
                                   aggregation=True, save_fig=False):
    
    target_file = data_extraction.file_selection(output_dir, temperature, deposition_rate)
    positions = data_extraction.read_variable_file(target_file, 'positions')
    
    plot_configuration(positions=positions, step_target=step_target, radius=radius, color_phi=False, diameter=diameter, 
                        x_limit=x_limit, y_limit=y_limit, aggregation=aggregation, 
                        save_fig=save_fig, save_name=f"config_T{temperature}_t{deposition_rate}.pdf",
                        set_title=f'Temperature={temperature}, Deposition Interval={deposition_rate}')

def animate_position_aggregation(output_dir, temperature, deposition_rate, line_length=0.4, area=50, color_p=False, frame_skip=1, frame_size=800, fps=20):
    target_file = data_extraction.file_selection(output_dir, temperature, deposition_rate)
    positions = data_extraction.read_variable_file(target_file, 'positions')
    
    animate_position_simple(positions=positions, output_name=f"animation_T{temperature}_t{deposition_rate}.mp4",
                            line_length=line_length, area=area, color_p=color_p, frame_skip=frame_skip, frame_size=frame_size, fps=fps)

############## Plot Heatmaps of parameters for output directory of aggregation with various T and dt ##############

def plot_parameter_heatmap(parameter_matrix, parameter_name, temperatures, deposition_rates, cmap='RdYlBu'):
    fig, ax = plt.subplots(figsize=(8, 6))
    
    c = ax.contourf(deposition_rates, temperatures, parameter_matrix, levels=100, cmap=cmap)

    fig.colorbar(c, ax=ax, label=parameter_name)
    ax.set_title(f'{parameter_name} Heatmap')
    ax.set_xlabel('Deposition Intervals')
    ax.set_ylabel('Temperature')

    plt.show()

