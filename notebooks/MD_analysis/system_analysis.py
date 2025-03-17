import numpy as np
import data_extraction

def calculate_radius_of_gyration(file_path):
    
    positions = data_extraction.read_variable_file(file_path, 'positions')['data']


    center_of_mass = np.mean(positions, axis=0)

    # Compute squared distances from center of mass
    squared_distances = np.sum((positions - center_of_mass) ** 2, axis=1)

    # Compute Radius of Gyration
    Rg = np.sqrt(np.mean(squared_distances))

    return Rg
