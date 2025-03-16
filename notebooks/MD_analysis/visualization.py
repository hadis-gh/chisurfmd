import numpy as np
import matplotlib.pyplot as plt

#  in the notebook you can try:
# %load_ext autoreload
# %autoreload 2

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

def plot_energies(kinetic_energy, potential_energy, m_temperatures, show_potential=True, temperature_show=False, temperature_label=False):
    kinetic_energy_data = kinetic_energy['data']
    temp_labels = kinetic_energy['temperature label']
    potential_energy_data = potential_energy['data']

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
        
def plot_order_parameter(order, m_temperatures, type, temperature_show=False, temperature_label=False):
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

    ax.set_ylabel("Order Parameter")
    ax.set_title(f"{type} Order Parameter")
        
    plt.tight_layout()
    plt.show()

    if temperature_show:
        plot_temperature(m_temperatures)

############## snapshot of system over steps/ single file --one MD simulation or one aggregation ##############

def plot_snapshot_configuration(positions, shot=-1, radius=0.4, color_phi=False, diameter=20, x_limit=(0, 50),y_limit=(0, 50), save_fig=False):

    x0, x1 = x_limit
    y0, y1 = y_limit
        
    positions_data = positions['data']
    availble_shot = positions_data.shape[0]
    print('available steps: ', availble_shot, '\tchoosed step: ', shot)

    x = positions_data[shot, :, 0]
    y = positions_data[shot, :, 1]
    if positions_data.shape[2] == 3:
        phi = positions_data[shot, :, 2]

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
            
            scatter = ax.scatter(
                x, y, 
                s=diameter,
                edgecolors='black',
                facecolors='none',
                alpha=0.8
            )
            x_end = x + radius * np.cos(phi)
            y_end = y + radius * np.sin(phi)
            
            for i in range(positions_data.shape[1]):
                ax.plot([x[i], x_end[i]], [y[i], y_end[i]], color='black', linewidth=0.8, alpha=0.8)

    ax.set_xlabel('X Position')
    ax.set_ylabel('Y Position')
    if shot == -1:
        ax.set_title(f"Particle's configuration - Last step")
    elif shot == 0:
        ax.set_title(f"Particle's configuration - First step")
    else:
        ax.set_title(f"Particle's configuration - {shot/availble_shot * 100:.0f}%")

    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)

    ax.set_aspect('equal', adjustable='box')
    ax.grid(linestyle='--', alpha=0.5)
    ax.set_axisbelow(True)

    plt.tight_layout()
    if save_fig:
        plt.savefig(f'config_{shot}.pdf')
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

        if step % (100 * frame_skip) == 0:
            print(f"Processing frame {step}/{num_steps}")

    out.release()
    print("Video saved as", output_name)

############## Depiction and alalysis of phi and delta phi of particles ##############

