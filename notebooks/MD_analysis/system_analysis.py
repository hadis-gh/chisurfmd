import numpy as np
import os
import re
import data_extraction
from scipy.spatial import ConvexHull

#------------------------------   POTENTIAL   ------------------------------

SIGMA = 1.0
EPSILON = 1.0

def lj_potential(r, sigma=SIGMA, epsilon=EPSILON):
    if np.any(r <= 0):
        raise ValueError("Distance r must be greater than zero.")
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

def lj_force(r, sigma=SIGMA, epsilon=EPSILON):
    if np.any(r <= 0):
        raise ValueError("Distance r must be greater than zero.")
    return 48 * epsilon * (sigma**12 / r**13 - 0.5 * sigma**6 / r**7)

def lj_mod_potential(r, delta_phi, phi_order, angular_scale, angular_order, sigma=SIGMA, epsilon=EPSILON, alpha=np.pi): 
    lj_term = lj_potential(r, sigma=SIGMA, epsilon=EPSILON)
    A = angular_scale / r**angular_order
    angular_term = A *(1 + np.cos(phi_order * delta_phi + alpha))
    
    return lj_term + angular_term

def lj_mod_force(r, delta_phi, phi_order, angular_scale, angular_order, sigma=SIGMA, epsilon=EPSILON):
    lj_term = lj_force(r, sigma=SIGMA, epsilon=EPSILON)
    
    A = angular_scale / r**angular_order
    dA_dr = -angular_order * A / r
    
    radial_force = lj_term - dA_dr *(1 + np.cos(phi_order * delta_phi))
    angular_force = phi_order * A * np.sin(phi_order * delta_phi)
    
    return radial_force, angular_force

def function_m(phi1, phi2, gamma, h1, h2, d1, d2, n, alpha):
    # delta_phi = h2*d2*(phi2 - gamma) - h1*d1*(phi1 - gamma)
    delta_phi = phi2 - phi1 - (h1*d1 - h2*d2) * gamma
    return np.cos(n * delta_phi + alpha)

def lj_chiral_potential(r, phi1, phi2, phi_order, gamma, potential_type, angular_scale, angular_order, sigma=SIGMA, epsilon=EPSILON, alpha=np.pi): 
    lj_term = lj_potential(r, sigma=sigma, epsilon=epsilon)
    A = angular_scale / r**angular_order
    
    if potential_type == 'RRUU':
        h1, h2, d1, d2 = 1, 1, 1, 1
    elif potential_type == 'RLUU':
        h1, h2, d1, d2 = 1, -1, 1, 1
    elif potential_type == 'RRUD':
        h1, h2, d1, d2 = 1, 1, 1, -1
    elif potential_type == 'RLUD':
        h1, h2, d1, d2 = 1, -1, 1, -1
    else:
        raise ValueError(f"Unknown potential type: {potential_type}")
    
    angular_term = A * (1 + function_m(phi1, phi2, gamma, h1, h2, d1, d2, phi_order, alpha))
    
    return lj_term + angular_term

#------------------------------   AGGREGATION   ------------------------------

################# Calculate Parameters for Cluster Analysis in Aggregation #################

def calculate_gyration_radius(positions_data):
    positions_xy = positions_data[:, :, :2]  # Ignore phi
    COM = np.mean(positions_xy, axis=1)
    Rg = []

    for step in range(positions_xy.shape[0]):
        squared_distance = np.sum((positions_xy[step] - COM[step])**2, axis=1)
        Rg.append(np.sqrt(np.mean(squared_distance)))
        
    return np.array(Rg)

def calculate_fractal_dimension(positions_data):
    positions_xy = positions_data[:, :, :2]
    fractal_dimensions = []

    for step in range(positions_xy.shape[0]):
        Df = 1
        fractal_dimensions.append(Df)

    return np.array(fractal_dimensions)

def calculate_convexHull_area(positions_data):
    positions_xy = positions_data[:, :, :2]
    perimeter_area = []

    for step in range(positions_xy.shape[0]):
        valid_particles = ~np.isnan(positions_xy[step]).any(axis=1)  # Ignore NaN values
        valid_positions = positions_xy[step][valid_particles]

        hull = ConvexHull(valid_positions)

        PA_ratio = hull.area 
        perimeter_area.append(PA_ratio)

    return np.array(perimeter_area)

def compute_parameter_matrix(output_dir, temperatures, deposition_rates, parameter_name):
    temperature_numbers = len(temperatures)
    deposition_rates_numbers = len(deposition_rates)
    
    parameter_matrix = np.zeros((temperature_numbers, deposition_rates_numbers))

    for i, temperature in enumerate(temperatures):
        for j, deposition_rate in enumerate(deposition_rates):
                        
            file_path = data_extraction.file_selection(output_dir, temperature, deposition_rate)
            
            positions_data = data_extraction.read_variable(file_path, 'positions', print_message=False)['data']
            
            if parameter_name == 'gyration radius':
                parameter_data = calculate_gyration_radius(positions_data)
            elif parameter_name == 'fractal dimention':
                parameter_data = calculate_fractal_dimension(positions_data)
            elif parameter_name == 'convexHull':
                parameter_data = calculate_convexHull_area(positions_data)
            parameter_value = parameter_data[-1]

            parameter_matrix[i, j] = parameter_value
            
            print(f"File {i+1}/{temperature_numbers}, {j+1}/{deposition_rates_numbers} stored!                 ", end='\r')

    return parameter_matrix

