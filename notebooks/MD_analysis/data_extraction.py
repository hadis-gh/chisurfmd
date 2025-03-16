import os
import numpy as np
from adios2 import Stream
from adios2 import FileReader



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

def read_variable_directory(output_dir, variable_name):
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
    
    print(f'read {variable_name} successfully!')        
    return {'data':np.array(variable_data), 'temperature label':np.array(temperatures)}

def read_variable_file(file_name, variable_name):
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
    
    print(f'read {variable_name} successfully!')        
    return {'data':np.array(variable_data), 'temperature label':np.array(temperatures)}

def read_variable(name, variable_name):
    if name.endswith(".bp"):
        return read_variable_file(name, variable_name)
    elif os.path.isdir(name):
        return read_variable_directory(name, variable_name)

################# read attributes from output file and directory #################

def read_attributes_file(file_name, attribute_name):
    
    with FileReader(file_name) as reader:
        if attribute_name in reader.available_attributes():
            attribute_data = reader.read_attribute(attribute_name)
        else:
            raise ValueError(f"Variable {attribute_name} not found in {file_name}")
    
    print(f'{attribute_name} is {attribute_data}!')        
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
    
