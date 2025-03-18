import numpy as np
import os
import re
import data_extraction

def calculate_gyration_radius(positions_data):
    positions_xy = positions_data[:, :, :2]  # Ignore phi
    COM = np.mean(positions_xy, axis=1)
    Rg = []

    for step in range(positions_xy.shape[0]):
        squared_distance = np.sum((positions_xy[step] - COM[step])**2, axis=1)
        Rg.append(np.sqrt(np.mean(squared_distance)))
        
    return np.array(Rg)

def compute_radius_matrix(output_dir, temperatures, deposition_rates):
    temperature_numbers = len(temperatures)
    deposition_rates_numbers = len(deposition_rates)
    
    radius_matrix = np.zeros((temperature_numbers, deposition_rates_numbers))

    for i, temperature in enumerate(temperatures):
        for j, deposition_rate in enumerate(deposition_rates):
            
            formatted_temp = f"{temperature:.2f}"
            formatted_dep_rate = f"{deposition_rate:.1f}"
                        
            file_path = os.path.join(output_dir, f'aggregated_{formatted_temp}_{formatted_dep_rate}.bp')
            
            positions_data = data_extraction.read_variable(file_path, 'positions', print_message=False)['data']
            radius_data = calculate_gyration_radius(positions_data)
            # radius_value = np.mean(radius_data)
            radius_value = radius_data[-1]

            radius_matrix[i, j] = radius_value
            print(f"File {i+1}/{temperature_numbers}, {j+1}/{deposition_rates_numbers} stored!                 ", end='\r')

    return radius_matrix

