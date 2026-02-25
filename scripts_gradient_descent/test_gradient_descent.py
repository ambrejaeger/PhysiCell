# Let's write a multilevel parameter estimation tool:
# PhysiCOOL: A generalized framework for model Calibration and Optimization Of modeLing projects
# David Hormuth,Ines Goncalves,Caleb Phillips, Sandhya Prabhakaran
# revised: 08/23/2021


# Import the necessary libraries
import subprocess
from pathlib import Path
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import physicool
import physicool_extension
import re
import os
import scripts_sensibility_analysis.membrane_apoptosis_functions as membrane_apoptosis_functions
from scripts_sensibility_analysis.sensibility_analysis_functions import *
from scripts_sensibility_analysis.saved_data_modifications import modify_cell_data
from scripts_sensibility_analysis.hill_response import read_value_in_settings_file

# what is the name of the compiled project?
PROJECT_NAME = './test_death'  
# Select the cell-specific variables you need to load in for optimization
VARIABLES = ['ID', 'position_x', 'position_y','position_z']
# where is the model output stored?
STORAGE_PATH = Path('output')
# where is the configuration file located?
CONFIG_PATH = Path('config/PhysiCell_settings.xml')
NUMBER_OF_CELLS = 11

def read_output(storage_path, variables):
    #This function will read in location and ID of each cell and return it as "cells_df"
    cells_through_time = []
    timesteps = get_timesteps(storage_path)
    for timestep in timesteps:
        # Read the data saved at each time point
        cells = physicool_extension.get_cell_data(timestep, storage_path, variables)
        number_of_cells = len(cells['ID'])

        # Store the data for each cell
        for i in range(number_of_cells):
            cells_data = [cells[variable][i] for variable in variables] + [timestep]
            cells_through_time.append(cells_data)
            
    cells_df = pd.DataFrame(cells_through_time, columns=['ID', 'x', 'y','z','time'])
    return cells_df


def get_timesteps(storage_path):
    """Returns the number of output XML files in the storage directory."""
    files = list(storage_path.glob('output*.xml'))
    number_of_output_files = len(files)
    if number_of_output_files == 0:
        return range(0)

    nums = []
    for f in files:
        m = re.search(r"output0*(\d+)\.xml$", f.name)
        if m:
            nums.append(int(m.group(1)))
        else:
            m2 = re.search(r"(\d+)", f.name)
            if m2:
                nums.append(int(m2.group(1)))

    if not nums:
        # fallback to simple 0..n-1
        return range(number_of_output_files)

    start = min(nums)
    return range(start, start + number_of_output_files)

def compute_cell_count(cells_df):
    # This function counts the total number of cells over time
    max_values = cells_df.max()
    min_values = cells_df.min()
    start_index = int(min_values['time'])
    last_index = int(max_values['time'])

    cell_count = np.zeros((last_index+1 - start_index))
    for n in range(start_index,last_index+1 ): # loops through time
        cell_count[n - start_index] = (cells_df['time'] == n).sum()      
    return cell_count


def run_simulation():
    command = PROJECT_NAME
    subprocess.run(command, shell=True) 
    
    
def run_pipeline(params):
    physicool_extension.update_config_file(params, CONFIG_PATH, 'epi_inter')
    run_simulation()
    cells = read_output(STORAGE_PATH, VARIABLES)
    avg_cc = compute_cell_count(cells).mean()
    label_file = os.path.join(STORAGE_PATH, 'initial.xml')
    mat_files = get_output_files(STORAGE_PATH)
    avg_pos_apop = membrane_apoptosis_functions.compute_median_apoptosis_position(label_file, mat_files)[1]
    apop_count = membrane_apoptosis_functions.compute_median_apoptosis_position(label_file, mat_files)[0]
    return avg_cc, avg_pos_apop

def run_pipeline_avg(params):
    # When calculating the objective function or the Jacobian, it may be beneficial to run the code...
    # several times and then average the replicates of multi-runs
    physicool_extension.update_config_file(params, CONFIG_PATH, 'epi_inter')
    run_simulation()
    cells = read_output(STORAGE_PATH, VARIABLES)
    avg_cc, lmh_quant, cell_count_rep_1 = compute_cell_count(cells)
    
    run_simulation()
    cells = read_output(STORAGE_PATH, VARIABLES)
    avg_cc_2, lmh_quant_2, cell_count_rep_2 = compute_cell_count(cells)
    
    run_simulation()
    cells = read_output(STORAGE_PATH, VARIABLES)
    avg_cc_3, lmh_quant_3, cell_count_rep_3 = compute_cell_count(cells)
    
    run_simulation()
    cells = read_output(STORAGE_PATH, VARIABLES)
    avg_cc_4, lmh_quant_4, cell_count_rep_4 = compute_cell_count(cells)


    avg_cc = (0.25)*(avg_cc+avg_cc_2+avg_cc_3+avg_cc_4);
    lmh_quant = (0.25)*(lmh_quant+lmh_quant_2+lmh_quant_3+lmh_quant_4)
    
    return avg_cc, lmh_quant

def update_saved_data(params, saved_data, config_file=CONFIG_PATH):
    path = "user_parameters/saving_folder"
    saved_data_folder = read_value_in_settings_file(config_file, path)
    if saved_data_folder:
        saved_cell_data_file = os.path.join(saved_data_folder, "cell_data.txt")
    else:
        raise ValueError("No saving folder found in the configuration file.")
    
    for key in params.keys():
        if key.split('/')[0] == 'cell':
            done = modify_cell_data(saved_cell_data_file, saved_data[key][0], saved_data[key][1], saved_data[key][2],params[key])
            if not done:
                raise ValueError(f"Modification of saved cell data failed for key {key} with value {params[key]}.")
    return        

def save_objective_function(objective_function, axis_0, axis_1, output_path, level_nbr=1):
    np.savez(f"{output_path}/objective_function_level_{level_nbr}.npz",axis_0=axis_0, axis_1=axis_1, objective_function=objective_function)

def read_objective_function(file_path):
    data = np.load(file_path)
    return data['axis_0'], data['axis_1'], data['objective_function']

if __name__ == "__main__":
    params = {"cell/model[@name='apoptosis']/death_rate" : 0.00025,
              "cell_rule/2/5" : 2.14,}
    save_path = {"cell/model[@name='apoptosis']/death_rate" :["type", 1, "Death/Model 0/Rate"]}

    storage_path = Path('output')
    variables = VARIABLES
    cells_df = read_output(storage_path, variables)
    label_file = os.path.join(STORAGE_PATH, 'initial.xml')
    mat_files = get_output_files(STORAGE_PATH)
    
    update_saved_data(params, save_path, config_file=CONFIG_PATH)

    """
    #This runs correctly
    key = "model[@name='apoptosis']/death_rate"
    #print(physicool_extension.get_cell_xml_stem(key, definition_name='epi_inter'))
    print(physicool_extension.update_config_file(params, CONFIG_PATH, 'epi_inter'))
    """
   
    
    # Parameter sweep 1 level
    number_of_levels = 4
    points_per_direction = 10
    percent_per_direction = 1
    num_params = len(params)
    about_point = np.array([0.0000532, 0.5])

    # Set lower and upperbound for model parameters
    param_lb = np.array([0.000001, 0.001]) # lower bound
    param_ub = np.array([0.0001, 1.0]) #upper bound 
    param_ub[0] = 0.008;

    parameters_in_sweep = np.zeros((num_params,1))
    if num_params == 1:
        objective_function = np.zeros((number_of_levels,points_per_direction))
        save_x = np.zeros((number_of_levels,points_per_direction))

    else:
        objective_function = np.zeros((number_of_levels,points_per_direction,points_per_direction))
        save_x = np.zeros((number_of_levels,points_per_direction))
        save_y = np.zeros((number_of_levels,points_per_direction))


    for n in range(number_of_levels):
        factor = percent_per_direction/(n*2+1)
        # Checks and make sure our parameters are within bounds, and generate parameter sweep
        if n == 0:
            param_1_sweep = np.linspace(param_lb[0],param_ub[0],points_per_direction)
        else:
            param_1_sweep = np.linspace(about_point[0]-factor*about_point[0],about_point[0]+factor*about_point[0],points_per_direction)
        
        param_1_sweep[param_1_sweep<param_lb[0]] = param_lb[0]
        param_1_sweep[param_1_sweep>param_ub[0]] = param_ub[0]
        save_x[n] = param_1_sweep

        if num_params>1:
            if n == 0:
                param_2_sweep = np.linspace(param_lb[1],param_ub[1],points_per_direction)
            else:
                param_2_sweep = np.linspace(about_point[1]-factor*about_point[1],about_point[1]+factor*about_point[1],points_per_direction)
            param_2_sweep[param_2_sweep<param_lb[1]] = param_lb[1]
            param_2_sweep[param_2_sweep>param_ub[1]] = param_ub[1]
            print(param_2_sweep)
            save_y[n] = param_2_sweep
            print(save_y[n])

        for a in range(points_per_direction):
            for b in range(points_per_direction):
                
                params = {"cell/model[@name='apoptosis']/death_rate" : param_1_sweep[a],
                "cell_rule/2/5" : param_2_sweep[b]}

                #update as well the cell data
                update_saved_data(params, save_path)

                avg_cc_model, pos_model = run_pipeline(params)
                print(avg_cc_model,pos_model)
                pos_data = 78.31404222264136 #Initial maximum position of cells along the y-axis
                cc_data = 726.0 #Initial number of cells

                if type(avg_cc_model) != type(cc_data) or type(pos_model) != type(pos_data):
                    objective_function[n][a][b] = None
                else:
                    objective_function[n][a][b] = ((avg_cc_model-cc_data)**2) + ((pos_model-pos_data)**2)
 
        I = np.argmin(objective_function[n])
        # I is optimal index, but in references to a points_per_direction X points_per_direcion
        x = int(np.floor(I/points_per_direction))
        #y = int(I-points_per_direction*x)
        
        
        about_point[0] = param_1_sweep[I]
        about_point[1] = param_2_sweep[I]
        print(about_point[0])
        print(about_point[1])

        output_path = "./grad_descent_output"
        save_objective_function(objective_function[n], param_1_sweep, param_2_sweep, output_path, level_nbr=1)