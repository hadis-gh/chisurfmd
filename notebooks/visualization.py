import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
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
    
    fig, axs = plt.subplots(figsize=(7, 7))
      
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
    
    fig, axs = plt.subplots(figsize=(7, 7))
      
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
    
    fig, axs = plt.subplots(figsize=(7, 7))
      
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
    
    fig, axs = plt.subplots(figsize=(7, 7))
      
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
    
#     fig, axs = plt.subplots(figsize=(7, 7))
      
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

def _apply_temp_xlabel(ax, n_steps, temp_labels, fontsize):
    ticks = np.linspace(0, n_steps, 10, dtype=int)
    tidx  = np.linspace(0, len(temp_labels) - 1, 10, dtype=int)
    ax.set_xlabel("Temperature", fontsize=fontsize)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"$T={temp_labels[i][0]:.1f}$" for i in tidx],
                       rotation=45, fontsize=fontsize - 2)

def plot_temperature(
    m_temperatures, temperature_label=False,
    colors=None, linewidth=1.5, fontsize=12,
    figsize=(5, 3), savepath=None,
):
    if colors is None:
        colors = _PUB_COLORS
    data = m_temperatures['data']
    temp_labels = m_temperatures['temperature label']

    fig, ax = plt.subplots(figsize=figsize)

    axis_labels = ['x', 'y', r'\omega']
    linestyles = ['-', '--', ':']

    for i, (lbl, ls) in enumerate(zip(axis_labels, linestyles)):
        ax.plot(data[:, :, i], label=f"$T_{{{lbl}}}$",
                color=colors[i], linewidth=linewidth)

    if temperature_label:
        _apply_temp_xlabel(ax, data.shape[0], temp_labels, fontsize)
    else:
        ax.set_xlabel("Steps", fontsize=fontsize)

    ax.set_ylabel("Measured Temperature", fontsize=fontsize)
    ax.set_title("Temperature Evolution", fontsize=fontsize + 1)
    ax.tick_params(axis='both', labelsize=fontsize - 2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.legend(bbox_to_anchor=(1, 1), loc='upper left',
              frameon=True, fancybox=False, shadow=False,
              fontsize=fontsize - 2, edgecolor='0.8')

    apply_style(ax, spine=False, grid=True, hide_top_right=False)
    plt.tight_layout()

    if savepath:
        fig.savefig(savepath, bbox_inches='tight')

    plt.show()
    return fig, ax
#------------------------------------------------------------
def plot_neighbors(
    neighbors_array, m_temperatures,
    temperature_show=False, temperature_label=False,
    colors=None, linewidth=1.5, fontsize=12,
    figsize=(5, 3), savepath=None,
):
    if colors is None:
        colors = ["#1F6E58", "#4B9D7B", "#90CDBC"]
    data = neighbors_array['data']
    temp_labels = neighbors_array['temperature label']

    fig, ax = plt.subplots(figsize=figsize)

    for i in range(data.shape[2]):
        ax.plot(data[:, :, i], label=f"Shell {i + 1}",
                color=colors[i % len(colors)],
                linewidth=linewidth)

    if temperature_label:
        _apply_temp_xlabel(ax, data.shape[0], temp_labels, fontsize)
    else:
        ax.set_xlabel("Steps", fontsize=fontsize)

    ax.set_ylabel("Number of Neighbors", fontsize=fontsize)
    ax.set_title("Neighbor Analysis", fontsize=fontsize + 1)
    ax.tick_params(axis='both', labelsize=fontsize - 2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.legend(bbox_to_anchor=(1, 1), loc='upper left',
              frameon=True, fancybox=False, shadow=False,
              fontsize=fontsize - 2, edgecolor='0.8')

    apply_style(ax, spine=False, grid=True, hide_top_right=False)
    plt.tight_layout()

    if savepath:
        fig.savefig(savepath, bbox_inches='tight')

    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

    return fig, ax
# ----------------------------------------------------------------------
def apply_style(ax, spine=False, grid=True, hide_top_right=True):
    """Apply common tick/grid/spine styling to an axis."""
    ax.tick_params(which='both', axis="both", direction="in")
    if grid:
        ax.grid(which='major', linestyle=":", alpha=.4, zorder=0)
        ax.grid(which='minor', linestyle=':', linewidth=0.4, alpha=.2, zorder=0)
        ax.set_axisbelow(True)
    if hide_top_right:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    if spine:
        for side in ['left', 'right', 'top', 'bottom']:
            ax.spines[side].set_linewidth(0.5)
            ax.spines[side].set_color('gray')

# ----------------------------------------------------------------------
_PUB_COLORS = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00']  # Wong colorblind-safe palette

def plot_energies(
    kinetic_energy, potential_energy, m_temperatures,
    show_potential=True, temperature_show=False, temperature_label=False,
    plot_window=None,
    colors=None,
    linewidth=1.5,
    fontsize=12,
    figsize=(5, 3),
    
    savepath=None,
):
    if colors is None:
        colors = _PUB_COLORS

    min_len = min(
        kinetic_energy['data'].shape[0],
        potential_energy['data'].shape[0],
        plot_window or int(1e9),
    )
    ke = kinetic_energy['data'][:min_len]
    pe = potential_energy['data'][:min_len]
    temp_labels = kinetic_energy['temperature label'][:min_len]

    fig, ax = plt.subplots(figsize=figsize)

    axis_labels = ['x', 'y', r'\omega']

    for i, lbl in enumerate(axis_labels):
        ax.plot(ke[:, :, i], label=f"$K_{{{lbl}}}$",
                color=colors[i], linewidth=linewidth)

    if show_potential:
        y2 = pe.flatten()
        ax.plot(y2, label='$U$',
                color=colors[3], linewidth=linewidth)
        ax.plot(y2 + np.sum(ke, axis=2)[:, 0], label='$E_\\mathrm{total}$',
                color=colors[4], linewidth=linewidth + 0.5)

    if temperature_label:
        n = ke.shape[0]
        ticks = np.linspace(0, n, 10, dtype=int)
        tidx  = np.linspace(0, len(temp_labels) - 1, 10, dtype=int)
        ax.set_xlabel("Temperature", fontsize=fontsize)
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"$T={temp_labels[i][0]:.1f}$" for i in tidx],
                           rotation=45, fontsize=fontsize - 2)
    else:
        ax.set_xlabel("Steps", fontsize=fontsize)

    ax.set_ylabel(r"Energy", fontsize=fontsize)
    ax.set_title("Energy Evolution", fontsize=fontsize + 1)
    ax.tick_params(axis='both', labelsize=fontsize - 2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.legend(bbox_to_anchor=(1, 1), loc='upper left',
              frameon=True, fancybox=False, shadow=False,
              fontsize=fontsize - 2, edgecolor='0.8')

    apply_style(ax, spine=False, grid=True, hide_top_right=False)
    plt.tight_layout()

    if savepath:
        fig.savefig(savepath, bbox_inches='tight')

    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

    return fig, ax

def plot_com_velocity(
    com_velocity, m_temperatures,
    temperature_show=False, temperature_label=False,
    colors=None, linewidth=1.5, fontsize=12,
    figsize=(5, 3), savepath=None,
):
    if colors is None:
        colors = _PUB_COLORS
    data = com_velocity['data']
    temp_labels = com_velocity['temperature label']

    fig, ax = plt.subplots(figsize=figsize)

    axis_labels = ['x', 'y']
    linestyles = ['-', '--']

    for i, (lbl, ls) in enumerate(zip(axis_labels, linestyles)):
        ax.plot(data[:, :, i], label=f"$v_{{{lbl}}}$",
                color=colors[i], linewidth=linewidth)

    if temperature_label:
        _apply_temp_xlabel(ax, data.shape[0], temp_labels, fontsize)
    else:
        ax.set_xlabel("Steps", fontsize=fontsize)

    ax.set_ylabel("CoM Velocity", fontsize=fontsize)
    ax.set_title("Center of Mass Velocity Evolution", fontsize=fontsize + 1)
    ax.tick_params(axis='both', labelsize=fontsize - 2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.legend(bbox_to_anchor=(1, 1), loc='upper left',
              frameon=True, fancybox=False, shadow=False,
              fontsize=fontsize - 2, edgecolor='0.8')

    apply_style(ax, spine=False, grid=True, hide_top_right=False)
    plt.tight_layout()

    if savepath:
        fig.savefig(savepath, bbox_inches='tight')

    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

    return fig, ax

def plot_com_ang_velocity(
    com_ang_velocity, m_temperatures,
    temperature_show=False, temperature_label=False,
    colors=None, linewidth=1.5, fontsize=12,
    figsize=(5, 3), savepath=None,
):
    if colors is None:
        colors = _PUB_COLORS
    data = com_ang_velocity['data']
    temp_labels = com_ang_velocity['temperature label']

    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(data.flatten(), color=colors[0], linewidth=linewidth)

    if temperature_label:
        _apply_temp_xlabel(ax, data.shape[0], temp_labels, fontsize)
    else:
        ax.set_xlabel("Steps", fontsize=fontsize)

    ax.set_ylabel("Angular Velocity", fontsize=fontsize)
    ax.set_title("Center of Mass Angular Velocity Evolution", fontsize=fontsize + 1)
    ax.tick_params(axis='both', labelsize=fontsize - 2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    apply_style(ax, spine=False, grid=True, hide_top_right=False)
    plt.tight_layout()

    if savepath:
        fig.savefig(savepath, bbox_inches='tight')

    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

    return fig, ax
        
def plot_order_parameter(
    order, m_temperatures,
    temperature_show=False, temperature_label=False,
    colors=None, linewidth=1.5, fontsize=12,
    figsize=(5, 3), savepath=None,
):
    if colors is None:
        colors = _PUB_COLORS
    data = order['data']
    temp_labels = order['temperature label']

    fig, ax = plt.subplots(figsize=figsize)

    ax.plot(data[:, :], color=colors[0], linewidth=linewidth)

    if temperature_label:
        _apply_temp_xlabel(ax, data.shape[0], temp_labels, fontsize)
    else:
        ax.set_xlabel("Steps", fontsize=fontsize)

    ax.set_ylabel("Orientation Order", fontsize=fontsize)
    ax.set_title("Orientation Order Parameter", fontsize=fontsize + 1)
    ax.tick_params(axis='both', labelsize=fontsize - 2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    apply_style(ax, spine=False, grid=True, hide_top_right=False)
    plt.tight_layout()

    if savepath:
        fig.savefig(savepath, bbox_inches='tight')

    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

    return fig, ax

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

    T = positions["temperature label"][0]
    
    fig, ax = plt.subplots(figsize=(5, 3.5))

    if positions_data.shape[2] == 3:
        phi = positions_data[step_target, :, 2]
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
    # else:
    #     ...
    if set_title:
        ax.set_title(set_title)
    else:
        if step_target == -1:
            ax.set_title(f"Last step")
        elif step_target == 0:
            ax.set_title(f"First step")
        else:
            ax.set_title(f"evolution: {step_target/availble_steps * 100:.0f}%")

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)

    ax.set_aspect('equal', adjustable='box')
    ax.set_axisbelow(True)
    apply_style(ax, spine=False, grid=True, hide_top_right=False)

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

def animate_position_simple(positions, output_name, line_length=0.4, area=50, color_p=False, frame_skip=1, frame_size=800, fps=20, patchNums=3):
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
                for i in range(patchNums):
                    # Draw direction line from center to end of radius
                    phi += 2 * np.pi /patchNums
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

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))

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

def plot_hist_phi_steps(positions, step_target, step_window=1000, bins_num=100):
    positions_data = positions['data']
    
    step_max = step_target + step_window
    print(f"available steps: {positions_data.shape[0]}, chosen step: {step_target} - {step_max}")
    
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize=(6, 4))

    ax.hist(positions_data[step_target:step_max,:,2].flatten(), bins=bins_num, color='#954965', alpha=0.7)

    ax.set_title(r'Histogram of particles $\phi$')
    ax.set_xlabel(r'Rotational Angle $\phi$')
    # ax.set_ylabel('Frequency')

    ax.set_axisbelow(True)
    ax.grid(color='gray', linestyle='dashed', alpha=0.5)

    ax.set_theta_zero_location("E")

    plt.tight_layout()
    plt.show()

def plot_deltaphi_hist_steps(positions, step_target, min_dis=20, step_window=1000, bins_num= 100):
    
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

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'},figsize=(6, 4))
    
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

def target_temperature_index (variable_name, target_temperature, cooling):
    data = variable_name['data']
    temperature_labels = variable_name['temperature label'].flatten()
    print(f"available temperatures are from {min(temperature_labels)} to {max(temperature_labels)}.")
          
    try:
        index = np.where(temperature_labels == target_temperature)[0][0 if cooling else 1]
    except IndexError:
        print(f"🞩🞩🞩 Target temperature {target_temperature} not found in temperature_labels. Please choose a valid temperature! 🞩🞩🞩")
        return 0

    total_snapshots = data.shape[0]
    return int(total_snapshots / len(temperature_labels) * index)

def plot_trajectory_temperature(positions, target_temperature, cooling=True, step_window=10000):
    positions_data = positions['data']

    step_min = target_temperature_index(positions, target_temperature, cooling)
    step_max = step_min + step_window

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))

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


def plot_snapshot_temperature(positions, handedness, 
                              target_temperature, cooling=True,
                              patchNums=None, line_length=0.45, 
                              area=20, radius=True, 
                              color_palette='hsv'):
    positions_data = positions['data']
    handedness_data = handedness['data']
    shot = target_temperature_index(positions, target_temperature, cooling)

    fig, ax = plt.subplots(figsize=(6, 4))
    
    x = positions_data[shot, :, 0]
    y = positions_data[shot, :, 1]
    phi = positions_data[shot, :, 2]  # orientation of each particle
    h = handedness_data[shot, :]

    # Define colors based on handedness
    colors = ['mistyrose' if val == 1 else 'lightsteelblue' for val in h]
    
    # Plot main particles
    ax.scatter(
        x, y, 
        s=110,
        edgecolors='black',
        facecolor=colors,
        alpha=0.8
    )

    if radius:
        if patchNums is not None:
            for a in range(patchNums):
                phi_a = phi + 2*np.pi * a / patchNums
                x_end = x + line_length * np.cos(phi_a)
                y_end = y + line_length * np.sin(phi_a)

                for i in range(len(x)):
                    ax.plot([x[i], x_end[i]], [y[i], y_end[i]], color='black', linewidth=0.8, alpha=0.8)
    # if radius:
    #     if patchNums is not None:
    #         phi = phi% (2*np.pi /patchNums)
    #         x_end = x + line_length * np.cos(phi)
    #         y_end = y + line_length * np.sin(phi)
    #         for i in range(len(x)):
    #             ax.plot([x[i], x_end[i]], [y[i], y_end[i]], color='black', linewidth=0.8, alpha=0.8)
                    
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

def plot_hist_phi_temperature(positions, target_t, step_window, cooling=True, bins_num=100):
    step_target = target_temperature_index(positions, target_t, cooling)
    plot_hist_phi_steps(positions, step_target, step_window, bins_num)

def plot_delta_phi_hist_temperature(positions, target_t, min_dis, step_window, cooling=True, bins_num=100):
    step_target = target_temperature_index(positions, target_t, cooling)
    plot_deltaphi_hist_steps(positions, step_target, min_dis, step_window, bins_num)

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

