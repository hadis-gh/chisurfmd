import numpy as np
import os
import re
import data_extraction
from scipy.spatial import ConvexHull

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

