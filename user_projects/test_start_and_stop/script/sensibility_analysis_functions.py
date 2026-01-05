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


def cell_type_height(cells_info, X_min, X_max, bin, cell_type):
    positions, cells_type = cells_info
    
    mask = cells_type == cell_type
    filtered_positions = positions[mask]
    bins = np.arange(X_min, X_max, bin)
    
    digitized = np.digitize(filtered_positions[:, 0], bins)
    heights = np.zeros((len(bins), 2))
    for i in range(len(bins)):
        mask = digitized == (i + 1)
        if np.any(mask):
            avg_y = np.mean(filtered_positions[mask, 1])
            heights[i] = [bins[i], avg_y]
    
    return heights


def cell_below(position, heights):
    mask = heights[:, 0] <= position[0]
    mask = np.argmax(heights[mask, 0])
    if heights[mask + 1, 0] >= position[0]:
        if heights[mask, 1] >= position[1]:
            return True
        else:
            return False
    else:
        print("This ain't right")
        return ValueError
    

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


def set_membrane_neighbors(neighbor_graph_file, cells_info, membrane_type):
    cells_type = cells_info[1]
    cells_ID = cells_info[2]

    mask = cells_type == membrane_type
    membrane_cells_ID = cells_ID[mask]

    neighbor_dict = {}

    with open(neighbor_graph_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                key, values = line.split(':', 1)
                if np.isin(int(key), membrane_cells_ID):
                    neighbors = np.asarray([int(x) for x in values.split(',')])
                    neighbor_dict[int(key)] = neighbors[np.isin(neighbors, membrane_cells_ID)]

    return neighbor_dict

def check_membrane_neighbors(cells_info, membrane_type, membrane_neighbors_0, X_min=-10000, X_max=10000):
    positions = cells_info[0]
    cells_type = cells_info[1]
    cells_ID = cells_info[2]

    mask = cells_type == membrane_type
    membrane_positions = positions[mask]
    x_filter_membrane = (membrane_positions[:,0] >= X_min) & (membrane_positions[:,0] <= X_max)
    membrane_positions = membrane_positions[x_filter_membrane]
    cells_ID = (cells_ID[mask])[x_filter_membrane]
    avg_dist = np.zeros(len(cells_ID))
    for i,cell_id in enumerate(cells_ID):
        #we are going to compute the distance between the cell and its neighbor at time t0
        #membrane_neighbors_0 is a dict with index the cell ID, and key a list of the initial membrane neighbors
        dist_neighbors = np.zeros(len(membrane_neighbors_0[cell_id]))
        for j,memb_neighbor in enumerate(membrane_neighbors_0[cell_id]):
            ind_neighbor = np.where(cells_ID == memb_neighbor)[0].item()

            dist_neighbors[j] = math.dist(membrane_positions[i], membrane_positions[ind_neighbor])

        avg_dist[i] = np.average(dist_neighbors)
    
    return avg_dist

def compute_membrane_integrity(mat_files, membrane_type, tolerance, X_min=-10000, X_max=10000):
    '''Compute for every saved file'''
    mat_file = "./temp_output/output00000001_cells.mat"
    initial_xml_file = "./temp_output/initial.xml"
    neighbor_graph_file = "./temp_output/output00000001_cell_neighbor_graph.txt"
    
    cells_info = extract_position_type_data(mat_file, initial_xml_file)
    membrane_neighbors_0 = set_membrane_neighbors(neighbor_graph_file, cells_info, membrane_type)
    avg_dist_1 = check_membrane_neighbors(cells_info, membrane_type, membrane_neighbors_0)

    count_breaks = [0.0 for _ in range(len(avg_dist_1))]
    for mat_file in mat_files[2:]:
        cells_info = extract_position_type_data(mat_file, initial_xml_file)
        avg_dist = check_membrane_neighbors(cells_info, membrane_type, membrane_neighbors_0)
        for i,dist in enumerate(avg_dist):
            if dist > avg_dist_1[i] + tolerance:
                
                count_breaks[i] = 1.0
    print(math.fsum(count_breaks))
    return count_breaks


def compute_membrane_integrity2(mat_files, temp_output_folder, membrane_type, tolerance, X_min=-10000, X_max=10000):
    '''Only compute at the last step of the run '''
    mat_file = os.path.join(temp_output_folder,"first_save_cells.mat")
    initial_xml_file = os.path.join(temp_output_folder,"initial.xml")
    neighbor_graph_file = os.path.join(temp_output_folder,"first_save_cell_neighbor_graph.txt")
    
    cells_info = extract_position_type_data(mat_file, initial_xml_file)
    membrane_neighbors_0 = set_membrane_neighbors(neighbor_graph_file, cells_info, membrane_type)
    avg_dist_1 = check_membrane_neighbors(cells_info, membrane_type, membrane_neighbors_0)

    count_breaks = [0.0 for _ in range(len(avg_dist_1))]
    cells_info = extract_position_type_data(mat_files[-1], initial_xml_file)
    avg_dist = check_membrane_neighbors(cells_info, membrane_type, membrane_neighbors_0)
    for i,dist in enumerate(avg_dist):
        if dist > avg_dist_1[i] + tolerance:
            count_breaks[i] = 1.0
    return count_breaks


def evaluate_membrane_integrity2(xml_file, param_treepaths, param_values, tolerance, save_output_folder, temp_output, restart=False, restart_int = 0, end_int = 0): 
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
            np.savetxt(os.path.join(output_folder,"param_values_membrane_integrity.txt"), param_values)  

        if os.path.exists(xml_file):
            shutil.copy(xml_file, temp_xml_file)
        
        start_file = restart_int
    else:
        if os.path.exists(xml_file) & (not os.path.exists(temp_xml_file)):
            shutil.copy(xml_file, temp_xml_file)

        param_values = np.loadtxt(os.path.join(output_folder,"param_values_membrane_integrity.txt"))
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
            modify_xml(temp_xml_file, param_treepaths[j][0], val, name_cell_def=param_treepaths[j][1], name_interact_cell_def=param_treepaths[j][2]) 

        #Now run the simulation that is loaded and made
        process0 = subprocess.run(
        ["./heterogeneity", temp_xml_file],
        capture_output=True,
        text=True,
        )

        #okay now we need to analyze the result
        mat_files = get_output_files(temp_output_folder)
        output_storage[i] = math.fsum(compute_membrane_integrity2(mat_files, temp_output_folder, 2, tolerance))
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


def restart_evaluate_membrane_integrity2(xml_file, tolerance, save_output_folder, temp_output):
    return evaluate_membrane_integrity2(xml_file, [], [], tolerance, save_output_folder, temp_output, restart=True)

def membrane_integrity(cells_info, membrane_type, membrane_thickness=0.0, X_min=-10000, X_max=10000):
    positions = cells_info[0]
    cells_type = cells_info[1]
    
    mask = cells_type == membrane_type
    membrane_positions = positions[mask]
    
    x_filter_membrane = (membrane_positions[:,0] >= X_min) & (membrane_positions[:,0] <= X_max)
    
    y_range_membrane = membrane_positions[x_filter_membrane,1].max() - membrane_positions[x_filter_membrane,1].min()
    print(y_range_membrane)
    if y_range_membrane > (membrane_thickness + 3):
        return y_range_membrane, False
    else:
        return y_range_membrane,True


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


def define_set_param(num_vars, names, bounds, groups=[], sample_size=32, seed=0): 
    print("Defining parameter space and sampling...")
    if len(groups) == 0:
        problem = {
            'num_vars': num_vars,
            'names': names,
            'bounds': bounds
        }
    else:
        problem = {
            'groups': groups,
            'num_vars': num_vars,
            'names': names,
            'bounds': bounds
        }

    param_values = sample(problem, sample_size, seed=seed) #Génère N*(2+D) jeux de paramètres avec D le nombre de paramètres et N un multiple de 2 fourni en argument
    return param_values

def evaluate_membrane_integrity(xml_file, param_treepaths, param_values): 
    
    temp_output_folder = os.path.join(os.getcwd(), "temp_output")
    temp_xml_file = os.path.join(temp_output_folder, os.path.basename(xml_file))

    if os.path.isdir(temp_output_folder):
        shutil.rmtree(temp_output_folder)
    
    os.makedirs(temp_output_folder, exist_ok=False)

    if os.path.exists(xml_file):
        shutil.copy(xml_file, temp_xml_file)

    output_storage = np.zeros(param_values.shape[0])
    memb_int = ["" for _ in range(param_values.shape[0])]
    #Changing output folder
    modify_xml(temp_xml_file, "save/folder", "temp_output")
    for i,X in enumerate(param_values): 
        #X is the set of values that should be modified in the xml
        for j, val in enumerate(X):
            modify_xml(temp_xml_file, param_treepaths[j][0], val, name_cell_def=param_treepaths[j][1], name_interact_cell_def=param_treepaths[j][2]) 
        
        #Now run the simulation that is loaded and made
        process0 = subprocess.run(
        ["./heterogeneity", temp_xml_file],
        capture_output=True,
        text=True,
        )

        #okay now we need to analyze the result
        initial_xml_file = os.path.join(temp_output_folder, "initial.xml")
        mat_files = get_output_files(temp_output_folder)
        initialized = False
        thickness = 0.0
        for mat_file in mat_files:
            cells_info = extract_position_type_data(mat_file, initial_xml_file)
            print(mat_file)
            if not initialized: 
                print("running")
                thickness, m = membrane_integrity(cells_info, 2)
                output_storage[i] = thickness
                initialized = True
            else:
                output_storage[i],m = membrane_integrity(cells_info, 2, membrane_thickness=thickness)
                if m == False:
                    memb_int[i] = mat_file
                    print(len(memb_int), " is the size of the list and the index is ", i)
        print("Run ", i, " is done")

        process1 = subprocess.run(
        ["make", "gif", "OUTPUT=temp_output"],
        capture_output=True,
        text=True
        )
        shutil.copyfile("./temp_output/out.gif", f"./output/out_{i}.gif")

    if os.path.isdir("./output"):
        with open('./output/param_names', 'w') as fp:
            for item in param_treepaths:
                fp.write("%s\n" % item)
        np.savetxt("./output/param_values_membrane_integrity.txt", param_values)
        np.savetxt("./output/membrane_integrity.txt", output_storage)
    return output_storage

def analyze_sobol(result_file, param_names_file, bounds):
    import ast
    with open(param_names_file, 'r') as f:
        names = [ast.literal_eval(line.strip()) for line in f if line.strip()]
    names = np.asarray([f"{';'.join(sublist)}" for sublist in names])
   
    print(type(names))
    print(type(bounds))

    problem = {
        'num_vars': len(param_names_file),
        'names': names,
        'bounds': bounds
    }
    
    with open(result_file, 'r') as f:
        lines = f.readlines()

    sorted_lines = sorted(lines, key=lambda x: int(x.split()[0]))
    Y = np.asarray([float(line.split()[1]) for line in sorted_lines])
    print("Array:", Y)
    print("Length:", len(Y))
    
    Si = analyze(problem, Y, print_to_console=True)
    return Si

def process_file(input_file, output_file):
    """
    Process the file to remove duplicate entries where first column is identical
    AND second column values are the same.
    """
    # Dictionary to store values: {first_col_value: [(second_col_value, original_line)]}
    data_dict = {}
    
    # Read and process the input file
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Parse all lines
    for line_num, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
            
        parts = line.split()
        if len(parts) >= 2:
            first_col = parts[0]
            second_col = parts[1]
            
            # Store in dictionary
            if first_col not in data_dict:
                data_dict[first_col] = []
            data_dict[first_col].append((second_col, line))
    
    # Prepare output lines
    output_lines = []
    
    # Process each unique first column value
    for first_col in sorted(data_dict.keys()):
        entries = data_dict[first_col]
        
        if len(entries) == 1:
            # No duplicates, keep the line
            output_lines.append(entries[0][1])
        else:
            # Check if all second column values are identical
            first_value = entries[0][0]
            all_same = all(entry[0] == first_value for entry in entries[1:])
            
            if all_same:
                # Keep only the first occurrence
                output_lines.append(entries[0][1])
                # Optional: print a message about duplicates removed
                print(f"Removed {len(entries)-1} duplicate(s) for {first_col} (all values = {first_value})")
            else:
                # Different values - keep all lines (or modify as needed)
                for entry in entries:
                    output_lines.append(entry[1])
                print(f"WARNING: Different values for {first_col}: {[e[0] for e in entries]}")
    
    # Write to output file
    with open(output_file, 'w') as f:
        for line in output_lines:
            f.write(line + '\n')
    
    # Count and display results
    original_count = sum(1 for line in lines if line.strip())
    new_count = len(output_lines)
    
    print(f"\nResults:")
    print(f"Original file: {original_count} lines")
    print(f"Processed file: {new_count} lines")
    print(f"Lines removed: {original_count - new_count}")
    print(f"Output saved to: {output_file}")
    
    return new_count

def main():
    
    param = [["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_basal", "epi_basal"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_basal", "epi_inter"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_inter", "epi_basal"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_basal", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "epi_basal"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "conjonctif"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "conjonctif", "membrane"]
    ]

    bounds = [[0.0, 1.0] * len(param)]
    names = ["_".join(p) for p in param]
    num_vars = len(param)

    param_values = define_set_param(num_vars, names, bounds, sample_size=128)
    print(param_values)
    print(len(param_values))
    print(os.getcwd())
    output_file = "./output_integrity3/membrane_integrity_processed.txt"
    

    param = [["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_basal", "epi_basal"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_basal", "epi_inter"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_basal", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "epi_basal"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "conjonctif"]
    ]
    
    result_file = "./output_integrity3/membrane_integrity.txt"

    process_file(result_file, output_file )
    param_names_file = "./output_integrity3/param_names.txt"
    bounds = [[0.0, 1.0] * len(param)]
    Si = analyze_sobol(result_file, param_names_file, bounds)
    Si.plot()
    plt.show()
    
if __name__ == '__main__':
    main()