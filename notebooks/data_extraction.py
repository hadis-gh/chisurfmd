import os
import re
import numpy as np
from adios2 import Stream
from adios2 import FileReader

#------------------------------   SINGLE MD   ------------------------------

################# print information of output file and directory #################

def print_output_file_info(file_path, print_summary=True):
    attributes = {}
    variables = {}

    with Stream(file_path, "r") as stream:
        for step_index, _ in enumerate(stream.steps()):
            if step_index == 0:
                for attr in stream.available_attributes():
                    attributes[attr] = stream.read_attribute(attr)

            for var in stream.available_variables():
                if var not in variables:
                    variables[var] = []
                variables[var].append(stream.read(var))

            if step_index > 10000:
                print("Warning: Too many steps, breaking the loop.")
                break    

    for var in variables:
        variables[var] = np.array(variables[var])

    if print_summary:
        print("Attributes Summary:")
        print("-------------------")
        for idx, (name, value) in enumerate(attributes.items(), start=1):
            print(f"{idx}. {name:<40} {value}")

        print("\nVariables Summary:")
        print("------------------")
        for idx, (name, arr) in enumerate(variables.items(), start=1):
            shape = arr.shape
            description = f"{shape[0]} * {list(shape[1:])}" if arr.ndim > 1 else shape[0]
            print(f"{idx}. {name:<35} {description}")

def print_output_dir_info(output_dir, print_summary=True, index=0):
    file_path = os.path.join(output_dir, f"run_{index}.bp")
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    print_output_file_info(file_path, print_summary)

def print_output_info(name):
    if name.endswith(".bp"):
        print_output_file_info(name)
    elif os.path.isdir(name):
        print_output_dir_info(name)    
    else:
        raise ValueError(f"Invalid input: {name}. It should be a .bp file or a directory.")

################# read variables from output file and directory #################

def read_variable_directory(output_dir, variable_name, print_message=True):
    pattern = re.compile("run_([0-9]+).bp")
    run_numbers = [int(pattern.match(x)[1]) for x in os.listdir(output_dir) if pattern.match(x)]
    
    if not run_numbers:
        raise ValueError(f"No run files found in {output_dir}")
    
    run_numbers.sort()
    variable_data = []
    temperatures = []
    
    for i, run_num in enumerate(run_numbers):
        file_path = os.path.join(output_dir, f"run_{run_num}.bp")
        with FileReader(file_path) as reader:
            if variable_name in reader.available_variables():
                var_info = reader.available_variables()[variable_name]
                steps = int(var_info.get("AvailableStepsCount", 1))
    
                data = reader.read(variable_name, step_selection=[0, steps])
                if (data.ndim)== 1:
                    data = np.reshape(data, (1*steps, data.shape[0]//steps))
                else:
                    data = np.reshape(data, (1*steps, data.shape[0]//steps, data.shape[1]))

                variable_data.extend(data)
                
                temp_label = reader.read_attribute("temperature")
                temp_label = temp_label.flatten()
                temperatures.append(temp_label)
                
            else:
                raise ValueError(f"Variable {variable_name} not found in {file_path}")
    if print_message:
        print(f'read {variable_name} successfully!', end='\r')        
    return {'data':np.array(variable_data), 'temperature label':np.array(temperatures)}

def read_variable_file(file_name, variable_name, print_message=True):
    variable_data = []
    temperatures = []
    
    with FileReader(file_name) as reader:
        if variable_name in reader.available_variables():
            var_info = reader.available_variables()[variable_name]
            steps = int(var_info.get("AvailableStepsCount", 1))

            data = reader.read(variable_name, step_selection=[0, steps])
            if (data.ndim)== 1:
                data = np.reshape(data, (1*steps, data.shape[0]//steps))
            else:
                data = np.reshape(data, (1*steps, data.shape[0]//steps, data.shape[1]))

            variable_data.extend(data)
            
            temp_label = reader.read_attribute("temperature")
            temp_label = temp_label.flatten()
            temperatures.append(temp_label)
            
        else:
            raise ValueError(f"Variable {variable_name} not found in {file_name}")
    if print_message:
        print(f'read {variable_name} successfully!', end='\r')        
    return {'data':np.array(variable_data), 'temperature label':np.array(temperatures)}

def read_variable(name, variable_name, print_message=True):
    if name.endswith(".bp"):
        return read_variable_file(name, variable_name, print_message)
    elif os.path.isdir(name):
        return read_variable_directory(name, variable_name, print_message)

################# read attributes from output file and directory #################

def read_attributes_file(file_name, attribute_name):
    
    with FileReader(file_name) as reader:
        if attribute_name in reader.available_attributes():
            attribute_data = reader.read_attribute(attribute_name)
        else:
            raise ValueError(f"Variable {attribute_name} not found in {file_name}")
    
    # print(f'{attribute_name} is {attribute_data}!')        
    return attribute_data

def read_attributes_directory(output_dir, attribute_name, index=0):
    file_path = os.path.join(output_dir, f"run_{index}.bp")
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    return read_attributes_file(file_path, attribute_name)

def read_attributes(name, attribute_name):
    if name.endswith(".bp"):
        return read_attributes_file(name, attribute_name)
    elif os.path.isdir(name):
        return read_attributes_directory(name, attribute_name)
    
################# read multiple variables from output file and directory #################

def read_variables_batch(name, variable_names, print_message=True):
    """
    Read multiple variables in a single pass through the data
    """
    if name.endswith(".bp"):
        return _read_variables_batch_file(name, variable_names, print_message)
    elif os.path.isdir(name):
        return _read_variables_batch_directory(name, variable_names, print_message)

def _read_variables_batch_file(file_name, variable_names, print_message=True):
    """
    Read multiple variables from a single file
    """
    result = {}
    
    with FileReader(file_name) as reader:
        # Read temperature attribute once
        temp_label = reader.read_attribute("temperature").flatten()
        
        for variable_name in variable_names:
            if variable_name in reader.available_variables():
                var_info = reader.available_variables()[variable_name]
                steps = int(var_info.get("AvailableStepsCount", 1))

                data = reader.read(variable_name, step_selection=[0, steps])
                if data.ndim == 1:
                    data = np.reshape(data, (1*steps, data.shape[0]//steps))
                else:
                    data = np.reshape(data, (1*steps, data.shape[0]//steps, data.shape[1]))

                result[variable_name] = {
                    'data': np.array(data),
                    'temperature label': temp_label
                }
            else:
                raise ValueError(f"Variable {variable_name} not found in {file_name}")
    
    if print_message:
        print(f'Read {len(variable_names)} variables successfully!')
    
    return result

def _read_variables_batch_directory(output_dir, variable_names, print_message=True):
    """
    Read multiple variables from a directory of files
    """
    pattern = re.compile("run_([0-9]+).bp")
    run_numbers = [int(pattern.match(x)[1]) for x in os.listdir(output_dir) if pattern.match(x)]
    
    if not run_numbers:
        raise ValueError(f"No run files found in {output_dir}")
    
    run_numbers.sort()
    
    # Initialize result structure
    result = {var: {'data': [], 'temperature label': []} for var in variable_names}
    
    for i, run_num in enumerate(run_numbers):
        file_path = os.path.join(output_dir, f"run_{run_num}.bp")
        with FileReader(file_path) as reader:
            # Read temperature attribute once per file
            temp_label = reader.read_attribute("temperature").flatten()
            
            for variable_name in variable_names:
                if variable_name in reader.available_variables():
                    var_info = reader.available_variables()[variable_name]
                    steps = int(var_info.get("AvailableStepsCount", 1))

                    data = reader.read(variable_name, step_selection=[0, steps])
                    if data.ndim == 1:
                        data = np.reshape(data, (1*steps, data.shape[0]//steps))
                    else:
                        data = np.reshape(data, (1*steps, data.shape[0]//steps, data.shape[1]))

                    result[variable_name]['data'].extend(data)
                    result[variable_name]['temperature label'].append(temp_label)
                else:
                    raise ValueError(f"Variable {variable_name} not found in {file_path}")
    
    # Convert lists to arrays
    for var in variable_names:
        result[var]['data'] = np.array(result[var]['data'])
        result[var]['temperature label'] = np.array(result[var]['temperature label'])
    
    if print_message:
        print(f'Read {len(variable_names)} variables from {len(run_numbers)} files successfully!')
    
    return result
#------------------------------   AGGREGATION   ------------------------------

################# MD aggregation directory,find file #################

def file_selection(output_dir, temperature, deposition_interval):
    formatted_temp = f"{temperature:.2f}"
    formatted_dep_rate = f"{deposition_interval:.1f}"
    return os.path.join(output_dir, f'aggregated_{formatted_temp}_{formatted_dep_rate}.bp')

def extract_ranges(output_dir):
    """Extracts temperature and deposition intervals values from filenames in the given directory."""
    
    pattern = re.compile(r"aggregated_(\d+\.\d+)_(\d+\.\d+)\.bp")
    temperatures = []
    deposition_rates = []
    
    for file in os.listdir(output_dir):
        match = pattern.match(file)
        if match:
            temperatures.append(float(match.group(1)))
            deposition_rates.append(float(match.group(2)))
    
    if not temperatures or not deposition_rates:
        print("No matching files found.")
        return {'temperatures': [], 'deposition_rates': []}

    temperatures = sorted(set(temperatures))
    deposition_rates = sorted(set(deposition_rates))

    print(f'Temperature range: \t{temperatures[0]} .. {temperatures[-1]}')
    print(f'Deposition Rate range: \t{deposition_rates[0]} .. {deposition_rates[-1]}')

    return {'temperatures': temperatures, 'deposition_rates': deposition_rates}

def extract_order_parameter_matrix(output_dir, temperatures, deposition_rates):
    
    temperature_numbers = len(temperatures)
    deposition_rates_numbers = len(deposition_rates)
    order_matrix = np.zeros((temperature_numbers, deposition_rates_numbers))

    for i, temperature in enumerate(temperatures):
        for j, deposition_rate in enumerate(deposition_rates):
            formatted_temp = f"{temperature:.2f}"
            formatted_dep_rate = f"{deposition_rate:.1f}"

            file_path = os.path.join(output_dir, f'aggregated_{formatted_temp}_{formatted_dep_rate}.bp')

            order_value = read_variable(file_path, 'orientational order')['data']

            order_matrix[i, j] = order_value[-1]

            print(f"File {i+1}/{temperature_numbers}, {j+1}/{deposition_rates_numbers} stored!                 ", end='\r')

    return order_matrix

def extract_neighbors_number_matrix(output_dir, temperatures, deposition_rates):
    
    temperature_numbers = len(temperatures)
    deposition_rates_numbers = len(deposition_rates)
    neighbors_matrix = np.zeros((temperature_numbers, deposition_rates_numbers))

    for i, temperature in enumerate(temperatures):
        for j, deposition_rate in enumerate(deposition_rates):
            formatted_temp = f"{temperature:.2f}"
            formatted_dep_rate = f"{deposition_rate:.1f}"

            file_path = os.path.join(output_dir, f'aggregated_{formatted_temp}_{formatted_dep_rate}.bp')

            neighbors = read_variable(file_path, 'number of neighbors')['data']

            neighbors_matrix[i, j] = neighbors[-1, :,-1]

            print(f"File {i+1}/{temperature_numbers}, {j+1}/{deposition_rates_numbers} stored!                 ", end='\r')

    return neighbors_matrix


