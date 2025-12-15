import os
import shutil
import glob
import scipy.io
import numpy as np
import xml.etree.ElementTree as ET
import re
import numpy as np
import matplotlib.pyplot as plt
from SALib.analyze.sobol import analyze
from SALib.sample.sobol import sample
import subprocess
from collections import defaultdict
import pandas as pd
import shutil
import math
import random

def parse_physicell_labels(xml_file):
    """Parse the XML file to get labels and metadata for each field, it works"""
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        # Find the labels section in the XML
        labels_elem = root.find(".//labels")
        if labels_elem is None:
            raise ValueError(f"Could not find labels in XML file: {xml_file}")
        
        # Parse each label
        labels = {}
        count = 0
        for label in labels_elem.findall('label'):
            index = int(label.get('index'))
            size = int(label.get('size'))
            units = label.get('units', 'none')
            name = label.text.strip() if label.text else f"field_{index}"
            
            # Sanitize name to ensure XML compatibility
            #name = ''.join(c if c.isalnum() or c in '_- ' else '_' for c in name)
            
            labels[name] = [np.arange(count, count+size), units]
            count += size
        return labels
    
    except Exception as e:
        print(f"Error parsing XML file {xml_file}: {str(e)}")
        return {}
    

def parse_physicell_metadata(xml_file):
    """Parse the XML file to get labels and metadata for each field, it works"""
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        # Find the labels section in the XML
        metadata = root.find(".//metadata")
        if metadata is None:
            raise ValueError(f"Could not find metadata in XML file: {xml_file}")
        
        # Parse each label
        metadata_dict = {}
        metadata_dict['current_time'] = metadata.find('current_time')
        metadata_dict['current_runtime'] = metadata.find('current_runtime')
           
        return metadata_dict
    
    except Exception as e:
        print(f"Error parsing XML file {xml_file}: {str(e)}")
        return {}


def modify_xml(xml_file, path, value, name_cell_def="", name_interact_cell_def=""):
    """Modify XML file with error handling for incorrect paths and names"""
    try:
        if not os.path.exists(xml_file):
            raise FileNotFoundError(f"XML file '{xml_file}' does not exist")
        
        modified = False
        tree = ET.parse(xml_file)
        root = tree.getroot()

        if name_cell_def == "":
            element = root.find(path)
            if element is None:
                raise ValueError(f"Path '{path}' not found in XML structure")
            element.text = str(value)
            modified = True

        else: 
            path_to_cell_def = path.split("cell_definition/", 2)[0]
            path_to_attribute = "./" + path.split("cell_definition/", 2)[1]
            
            if root.find(path_to_cell_def) is None:
                raise ValueError(f"Path: '{path_to_cell_def}' not found in XML structure")
            
            subroot = None
            for cell_def in root.findall(os.path.join(path_to_cell_def,"cell_definition")):
                #print(cell_def.get("name"))
                if cell_def.get("name") == name_cell_def:
                    subroot = cell_def
                    break
            if subroot is None:
                raise ValueError(f"Cell definition: '{name_cell_def}' not found in XML structure")
            
            if name_interact_cell_def == "":
                element = subroot.find(path_to_attribute)
                if element is None:
                    raise ValueError(f"Attribute path '{path_to_attribute}' not found in cell definition '{name_cell_def}'")
                element.text = str(value)
                modified = True
            else:

                target_interaction = None
                for interact_cell_def in subroot.findall(path_to_attribute):
                    if interact_cell_def.get("name") == name_interact_cell_def:
                        target_interaction = interact_cell_def
                        break
                
                if target_interaction is None:
                    error_msg = f"Interaction cell definition '{name_interact_cell_def}' not found in cell definition '{name_cell_def}' at path '{path_to_attribute}'"
                    raise ValueError(error_msg)
                
                target_interaction.text = str(value)
                modified = True

        if modified:
            tree.write(xml_file)
            return True
        else:
            return False

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return False
    except ET.ParseError as e:
        print(f"Error: Invalid XML format in file '{xml_file}': {e}")
        return False
    except ValueError as e:
        print(f"Error: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error: {e}")
        return False


def extract_position_type_data(mat_file, xml_file):

    try:
        labels = parse_physicell_labels(xml_file)
        if not labels:
            error_msg = f"Impossible to parse_physicell_labels of file at path {xml_file}"
            raise ValueError(error_msg)
        
        # Load mat file
        data = scipy.io.loadmat(mat_file)
        cells_data = data['cells']
        
        print(f"Processing {os.path.basename(mat_file)}")
        print(f"Cells array shape: {cells_data.shape}")
        
        num_cells = cells_data.shape[1]

        positions = np.empty((num_cells,3))
        cells_type = np.empty(num_cells)
        cells_ID = np.empty(num_cells)

        position_indices = labels["position"][0] 
        id_index = labels["ID"][0]
        cell_type_index = labels["cell_type"][0]
        
        for i in range(num_cells):
            try:
                positions[i] = cells_data[position_indices,i]
                cells_type[i] = cells_data[cell_type_index,i].item()
                cells_ID[i] = cells_data[id_index, i].item()
                
            except (IndexError, ValueError) as e:
                print(f"Warning: Error processing cell {i}: {e}")
                continue

        return positions, cells_type, cells_ID
        
    except Exception as e:
        print(f"Error processing file {os.path.basename(mat_file)}: {str(e)}")


def track_y_pos(cells_info):
    y_pos = np.zeros(cells_info.size())
    for cell in cells_info:
        if cell[1] == 1:
            y_pos[int(cell[2])] = cell[0][1]
    return y_pos


def cell_type_height(cells_info, X_min, X_max, bin_width, cell_type):
    positions = cells_info[0]
    cells_type = cells_info[1]
    mask = cells_type == cell_type
    filtered_positions = positions[mask]

    bins = np.arange(X_min, X_max, bin_width)
    digitized = np.digitize(filtered_positions[:, 0], bins)
    heights = np.zeros((len(bins), 2))
    
    for i in range(len(bins)):
        mask = digitized == (i + 1)
        
        if np.any(mask):
            avg_y = np.mean(filtered_positions[mask, 1])
            heights[i] = [bins[i], avg_y]
        else:
            heights[i] = [bins[i], np.nan]  # or 0
    return heights


def cell_below(position, heights):
    mask = heights[:, 0] < position[0]
    mask = np.argmax(heights[mask, 0])
    if mask + 1 < len(heights):
        if heights[mask + 1, 0] >= position[0]:
            if heights[mask, 1] >= position[1]:
                return True
            else:
                return False
         else:
        print("This ain't right")
        return ValueError
    else:
        return False
   
    

def cell_above(position, heights):
    mask = heights[:, 0] <= position[0]
    mask = np.argmax(heights[mask, 0])
    if heights[mask + 1, 0] >= position[0]:
        if heights[mask, 1] <= position[1]:
            return True
        else:
            return False
    else:
        print("This ain't right")
        return ValueError


def compute_membrane_permeability(mat_files, temp_output_folder, membrane_type, tolerance, X_min=-300, X_max=300):
    initial_xml_file = os.path.join(temp_output_folder,"initial.xml")
    for mat_file in mat_files: 
        cells_info = extract_position_type_data(mat_file, initial_xml_file)
        cells_type = cells_info[1] 
        positions = cells_info[0] 

        heights = cell_type_height(cells_info, X_min, X_max, 15, membrane_type)
        attracted = cells_type == 0
        result = cell_below(positions[attracted][0], heights)
        if result:
            print(mat_file)
            return 1
        
    return 0

def evaluate_membrane_permeability(xml_file, param_treepaths, param_values, tolerance, save_output_folder, temp_output, restart=False, restart_int = 0, end_int = 0): 
    output_folder = os.path.join(os.getcwd(), save_output_folder)
    save_output = os.path.join(output_folder, "run_output.txt")
    output_storage_file = os.path.join(output_folder,"membrane_integrity.txt")
    temp_output_folder = os.path.join(os.getcwd(), temp_output)
    temp_xml_file = os.path.join(temp_output_folder, os.path.basename(xml_file))
    start_file = 0

    if not os.path.isdir(temp_output_folder):
        os.makedirs(temp_output_folder, exist_ok=False)
    
    if end_int == 0:
        end_file = param_values.shape[0]
    else:
        end_file = end_int

    if not restart:
        if not os.path.isdir(output_folder):
            os.makedirs(output_folder, exist_ok=False)

        if os.path.isdir(output_folder): 
            with open(os.path.join(output_folder, "param_names.txt"), 'w') as fp:
                for item in param_treepaths:
                    fp.write("%s\n" % item)
            np.savetxt(os.path.join(output_folder,"param_values_membrane_permeability.txt"), param_values)  

        if os.path.exists(xml_file):
            shutil.copy(xml_file, temp_xml_file)
        
        start_file = restart_int
    else:
        param_values = np.loadtxt(os.path.join(output_folder,"param_values_membrane_permeability.txt"))
        pattern = r'["\'](.*?)["\']'
        with open(os.path.join(os.path.join(output_folder,"param_names.txt")), "r") as file:
            param_treepaths = [re.findall(pattern, line.strip()) for line in file if line.strip()]
        if restart_int == 0:
            #get_output_files returns a list of files with corresponding prefix and suffix in folder sorted by the run number in their name
            files = get_output_files(output_folder, prefix="out_", suffix=".gif")
            last_file = os.path.basename(files[-1])
            last_file_int = int(last_file[4: len(last_file) - 4])
            start_file = last_file_int + 1
        else:
            start_file = restart_int
    if end_file < start_file:
        print(f"The end index({end_file}) is inferior to the start index,({start_file}) check your input values")
        raise ValueError
            
    if not os.path.isfile(save_output):
        with open(save_output, 'w') as f:
                pass
    if not os.path.isfile(output_storage_file):
        with open(output_storage_file, 'w') as f:
            pass

    output_storage = np.zeros(param_values.shape[0])
    #Changing output folder
    modify_xml(temp_xml_file, "save/folder", temp_output)
    for i in range(start_file, end_file): 
        for j, val in enumerate(param_values[i,:]):
            if len(param_treepaths[j]) == 3:
                modify_xml(temp_xml_file, param_treepaths[j][0], val, name_cell_def=param_treepaths[j][1], name_interact_cell_def=param_treepaths[j][2]) 
            elif len(param_treepaths[j]) == 2:
                modify_xml(temp_xml_file, param_treepaths[j][0], val, name_cell_def=param_treepaths[j][1])
            elif len(param_treepaths[j]) == 1:
                modify_xml(temp_xml_file, param_treepaths[j][0], val)
            else:
                raise ValueError

       #Now run the simulation that is loaded and made
        process0 = subprocess.run(
        ["./heterogeneity", temp_xml_file],
        capture_output=True,
        text=True,
        )

        #okay now we need to analyze the result
        mat_files = get_output_files(temp_output_folder)
        output_storage[i] = compute_membrane_permeability(mat_files, temp_output_folder, 2, tolerance)
        print("Run ", i, " is done")

        process1 = subprocess.run(
        ["make", "gif", f"OUTPUT={temp_output_folder}"],
        capture_output=True,
        text=True
        )
        
        shutil.copyfile(f"{temp_output_folder}/out.gif", f"{output_folder}/out_{i}.gif")
        with open(output_storage_file, "a") as f:
            f.write(f"{i} {output_storage[i]}\n")
        with open(save_output, "a") as f:
            f.write(process0.stdout)
        with open(save_output, "a") as f:
            f.write(process1.stdout)

    #delete temp output at the end of the run
    if os.path.isdir(temp_output_folder):
            shutil.rmtree(temp_output_folder)

    return output_storage


def plot_cells_2D(positions):
    mask_x0 = positions[:, 0] != 0.0
    x = positions[mask_x0,0]
    
    y = positions[mask_x0, 1]
    plt.figure(figsize=(12, 6))
    print("y:", y)
    # Plot line
    plt.plot(x, y, 'b-', linewidth=3, label='data')

    # Customize the plot
    plt.xlabel('X values', fontsize=12)
    plt.ylabel('Y values', fontsize=12)
    plt.grid(True, alpha=0.3, linestyle='--')
    plt.legend(fontsize=10)
    plt.tight_layout()

    # Show the plot
    plt.show()


def get_output_files(path, prefix='output', suffix='_cells.mat'):
    
    pattern = os.path.join(path, prefix + "*" + suffix)
    files = glob.glob(pattern)

    pattern_regex = re.compile(rf'{re.escape(prefix)}(\d+){re.escape(suffix)}$')
    
    def extract_number(filepath):
        basename = os.path.basename(filepath)
        match = pattern_regex.search(basename)
        if match:
            return int(match.group(1))
        return -1 
    
    files.sort(key=extract_number)
    
    return files

def concatenate_results(files_mat, files_xml):
    return


def define_set_param(num_vars, names, bounds, seed=0): 
    print("Defining parameter space and sampling...")
    problem = {
        'num_vars': num_vars,
        'names': names,
        'bounds': bounds
    }

    param_values = sample(problem, 32,seed=seed) #Génère N*(2+D) jeux de paramètres avec D le nombre de paramètres et N un multiple de 2 fourni en argument
    return param_values


def analyze_sobol(xml_file, param_treepaths, param_values, num_vars, names, bounds):
    problem, param_values = define_set_param(num_vars, names, bounds)
    Y = evaluate_membrane_integrity(xml_file, param_treepaths, param_values)
    Si = analyze(problem, Y, print_to_console=True)
    return Si

def main():
    param = [["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "attracted"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "attracted", "membrane"],
    ["cell_definitions/cell_definition/phenotype/motility/speed", "attracted"],
    ["cell_definitions/cell_definition/phenotype/motility/migration_bias", "attracted"],
    ["cell_definitions/cell_definition/phenotype/volume/total", "attracted"]
    ]

    bounds = [[0.0, 1.0],
          [0.0, 1.0],
          [0.0, 1.0],
          [0.1, 5.0],
          [0.0, 1.0],
          [3000, 5500]]
    names = ["_".join(p) for p in param]
    num_vars = len(param)
    xml_file = "./config/PhysiCell_settings.xml"
    param_values = define_set_param(num_vars, names, bounds)

    param_values2 = define_set_param(num_vars, names, bounds)

    print(param_values[0])
    print(param_values2[0])

if __name__ == "__main__":
    main()
