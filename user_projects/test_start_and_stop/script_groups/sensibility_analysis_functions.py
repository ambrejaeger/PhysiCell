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

def define_problem(num_vars, names, bounds, groups=[]):
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
    return problem

def define_set_param(num_vars, names, bounds, groups=[], sample_size=32, seed=0): 
    print("Defining parameter space and sampling...")
    problem = define_problem(num_vars, names, bounds, groups)

    param_values = sample(problem, sample_size, seed=seed) #Génère N*(2+D) jeux de paramètres avec D le nombre de paramètres et N un multiple de 2 fourni en argument
    return param_values

def analyze_sobol(result_file, param_names_file, bounds, groups=[]):
    import ast
    with open(param_names_file, 'r') as f:
        names = [ast.literal_eval(line.strip()) for line in f if line.strip()]
    names = np.asarray([f"{';'.join(sublist)}" for sublist in names])

    problem = define_problem(len(names), names, bounds, groups)
    print(problem)
    
    with open(result_file, 'r') as f:
        lines = f.readlines()

    sorted_lines = sorted(lines, key=lambda x: int(x.split()[0]))
    Y = np.asarray([float(line.split()[1]) for line in sorted_lines])
    print("Array:", Y)
    print("Length:", len(Y))
    
    Si = analyze(problem, Y, print_to_console=True)
    return Si

def combine_files_with_header(file1_path, file2_path, output_path, column_names):
    
    df1 = pd.read_csv(file1_path, sep=r'\s+', header=None, names=['row_idx', 'nbr_breaks'])
    df2 = pd.read_csv(file2_path, sep=r'\s+', header=None)
        
    if len(column_names) == df2.shape[1]:
        df2.columns = column_names
    
    else:
        print(f"Warning: column_names has {len(column_names)} items, but file2 has {df2.shape[1]} columns")
        df2.columns = [f"col_{i+1}" for i in range(df2.shape[1])]
        
    nbr_breaks_list = [0.0] * len(df2)
        
    for idx, row in df1.iterrows():
        row_idx = int(row['row_idx'])  
        if 0 <= row_idx < len(df2):  
            nbr_breaks_list[row_idx] = row['nbr_breaks']
        
    df_combined = pd.DataFrame({'nbr_breaks': nbr_breaks_list})
    df_combined = pd.concat([df_combined, df2.reset_index(drop=True)], axis=1)
        
    df_combined.to_csv(output_path, sep='\t', index=False, float_format='%.6e')
        
    print(f"Combined file created: {output_path}")
    print(f"Shape of combined data: {df_combined.shape}")
    print(f"Rows from file1 assigned: {len(df1[df1['row_idx'] < len(df2)])}")
        
    return df_combined


def extract_X_Y(file_path, X_columns, Y_column='nbr_breaks'):

    df = pd.read_csv(file_path, sep='\t')
    if Y_column not in df.columns:
        raise ValueError(f"Y column '{Y_column}' not found in file. Available columns: {list(df.columns)}")
    
    Y = df[Y_column].values
    if isinstance(X_columns, str):
        if X_columns not in df.columns:
            raise ValueError(f"Column '{X_columns}' not found.")
        X = df[X_columns].values.tolist()
        
    elif isinstance(X_columns, int):
        if X_columns >= len(df.columns) or X_columns < 0:
            raise ValueError(f"Column index {X_columns} out of range.")
        X = df.iloc[:, X_columns].values.tolist()
        
    elif isinstance(X_columns, list):
        if not X_columns:
            raise ValueError("X_columns list is empty.")
        
        if all(isinstance(x, str) for x in X_columns):
            missing_cols = [col for col in X_columns if col not in df.columns]
            if missing_cols:
                raise ValueError(f"Columns not found: {missing_cols}")

            if len(X_columns) == 1:
                X = df[X_columns[0]].values.tolist()
            else:
                X = df[X_columns].values.T.tolist()
        
        elif all(isinstance(x, int) for x in X_columns):
            valid_indices = [idx for idx in X_columns if 0 <= idx < len(df.columns)]
            if len(valid_indices) != len(X_columns):
                raise ValueError(f"Some indices are out of range. Valid range: 0-{len(df.columns)-1}")
            
            if len(X_columns) == 1:
                X = df.iloc[:, X_columns[0]].values.tolist()
            else:
                X = df.iloc[:, X_columns].values.T.tolist()
        
        else:
            raise ValueError("X_columns list must contain all strings or all integers.")
    
    else:
        raise TypeError("X_columns must be string, integer, or list of strings/integers.")
    
    print(f"Successfully extracted:")
    print(f"  Y shape: {len(Y)} values")
    print(f"  X shape: {len(X) if isinstance(X[0], list) else '1D'} features")
    
    return X, Y



def plot_scatter_sets(Y, X, set_names=None, title="Scatter Plot", 
                     xlabel="Features", ylabel="nbr_breaks", 
                     figsize=(10, 6), save_path=None, show_plot=True):
    
    if isinstance(X, list):
        if X and isinstance(X[0], list):
            num_sets = len(X)
            for i, x_set in enumerate(X):
                if len(x_set) != len(Y):
                    raise ValueError(f"Set {i} has {len(x_set)} values, but Y has {len(Y)} values")
        else:
            num_sets = 1
            X = [X] 
    else:
        raise TypeError("X must be a list or list of lists")
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Default set names if not provided
    if set_names is None:
        set_names = [f"Feature Set {i+1}" for i in range(num_sets)]
    elif len(set_names) != num_sets:
        print(f"Warning: {len(set_names)} names provided for {num_sets} sets. Using default names.")
        set_names = [f"Feature Set {i+1}" for i in range(num_sets)]
    
    for i in range(num_sets):
        x_values = X[i]
        ax.scatter(x_values, Y,  
                  alpha=0.6, 
                  edgecolors='w', 
                  linewidth=0.5,
                  label=set_names[i])

    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {save_path}")
    
    if show_plot:
        plt.show()
    
    return fig, ax



import scipy.stats as stats

def plot_histogram_with_analysis(data, bins=20, title="Histogram with Analysis", 
                                figsize=(12, 5), show_plot=True):
    
    # Convert to numpy array
    data = np.asarray(data).flatten()
    
    # Remove any non-positive values for exponential comparison
    # (Exponential is only defined for positive values)
    data_positive = data[data > 0]
    
    # Calculate basic statistics
    stats_dict = {
        'n': len(data),
        'n_positive': len(data_positive),
        'mean': np.mean(data),
        'std': np.std(data),
        'median': np.median(data),
        'min': np.min(data),
        'max': np.max(data)
    }
    
    # Calculate Coefficient of Variation (CV)
    # CV = standard deviation / mean
    if stats_dict['mean'] != 0:
        cv = stats_dict['std'] / abs(stats_dict['mean'])
        stats_dict['cv'] = cv
        # Rule of thumb: CV > 1 suggests heavy-tailedness
        stats_dict['cv_heavy_tailed'] = cv > 1.0
    else:
        stats_dict['cv'] = float('inf')
        stats_dict['cv_heavy_tailed'] = True
    
    # Create figure with subplots
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # --- Plot 1: Histogram ---
    ax1 = axes[0]
    ax1.hist(data, bins=bins, color='skyblue', edgecolor='black', alpha=0.8, density=True)
    
    # Add vertical lines for mean and median
    ax1.axvline(stats_dict['mean'], color='red', linestyle='--', 
                linewidth=2, label=f'Mean = {stats_dict["mean"]:.3f}')
    ax1.axvline(stats_dict['median'], color='green', linestyle=':', 
                linewidth=2, label=f'Median = {stats_dict["median"]:.3f}')
    
    ax1.set_xlabel('Value', fontsize=12)
    ax1.set_ylabel('Density', fontsize=12)
    ax1.set_title('Histogram with Mean/Median', fontsize=13)
    ax1.legend(loc='upper right')
    ax1.grid(True, alpha=0.3)
    
    # Add statistics text to histogram
    stats_text = f"n = {stats_dict['n']}\n"
    stats_text += f"Mean = {stats_dict['mean']:.3f}\n"
    stats_text += f"Std = {stats_dict['std']:.3f}\n"
    stats_text += f"CV = {stats_dict['cv']:.3f}\n"
    stats_text += f"CV > 1: {stats_dict['cv_heavy_tailed']}"
    
    ax1.text(0.02, 0.98, stats_text,
             transform=ax1.transAxes,
             fontsize=10,
             verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # --- Plot 2: QQ-Plot against Exponential ---
    ax2 = axes[1]
    
    if len(data_positive) >= 10:
        # Sort the positive data
        sorted_data = np.sort(data_positive)
        n = len(sorted_data)
        
        # Calculate theoretical exponential quantiles
        # Exponential CDF: F(x) = 1 - exp(-λx)
        # Quantile function: Q(p) = -log(1-p)/λ
        # For standard exponential (λ=1): Q(p) = -log(1-p)
        
        # Empirical probabilities (avoiding 0 and 1)
        p = (np.arange(1, n + 1) - 0.5) / n
        theoretical_quantiles = -np.log(1 - p)  # Standard exponential (λ=1)
        
        # Scale theoretical quantiles to match data scale
        # Using method of moments: λ_hat = 1/mean
        lambda_hat = 1 / np.mean(data_positive)
        theoretical_quantiles = theoretical_quantiles / lambda_hat
        
        # Plot QQ-plot
        ax2.scatter(theoretical_quantiles, sorted_data, alpha=0.6, s=30, color='darkorange')
        
        # Add reference line (y = x)
        max_val = max(np.max(theoretical_quantiles), np.max(sorted_data))
        ax2.plot([0, max_val], [0, max_val], 'r--', linewidth=2, 
                label='Exponential reference', alpha=0.7)
        
        # Calculate correlation for fit assessment
        correlation = np.corrcoef(theoretical_quantiles, sorted_data)[0, 1]
        
        ax2.set_xlabel('Theoretical Exponential Quantiles', fontsize=12)
        ax2.set_ylabel('Sample Quantiles', fontsize=12)
        ax2.set_title(f'QQ-Plot vs Exponential (r = {correlation:.3f})', fontsize=13)
        ax2.legend(loc='upper left')
        ax2.grid(True, alpha=0.3)
        
        # Add interpretation text
        if correlation > 0.98:
            fit_quality = "Excellent exponential fit"
        elif correlation > 0.95:
            fit_quality = "Good exponential fit"
        elif correlation > 0.90:
            fit_quality = "Moderate exponential fit"
        else:
            fit_quality = "Poor exponential fit"
        
        # Check for heavy tails (points above line in right tail)
        # Use last 20% of points
        tail_start = int(0.8 * n)
        if tail_start < n - 1:
            tail_theoretical = theoretical_quantiles[tail_start:]
            tail_actual = sorted_data[tail_start:]
            tail_ratio = np.mean(tail_actual / tail_theoretical)
            
            if tail_ratio > 1.1:
                tail_info = "Heavier tails than exponential"
            elif tail_ratio < 0.9:
                tail_info = "Lighter tails than exponential"
            else:
                tail_info = "Similar tails to exponential"
            
            ax2.text(0.02, 0.98, f"{fit_quality}\n{tail_info}",
                     transform=ax2.transAxes,
                     fontsize=9,
                     verticalalignment='top',
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        
        stats_dict['qq_correlation'] = correlation
        stats_dict['data_positive_mean'] = np.mean(data_positive)
        stats_dict['lambda_hat'] = lambda_hat
    
    else:
        # Not enough positive data for QQ-plot
        ax2.text(0.5, 0.5, f"Not enough positive data\nfor QQ-plot\n(n_positive = {len(data_positive)})",
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax2.transAxes,
                 fontsize=12)
        ax2.set_title('QQ-Plot vs Exponential', fontsize=13)
        ax2.grid(True, alpha=0.3)
        stats_dict['qq_correlation'] = None
    
    # Overall title
    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    if show_plot:
        plt.show()
    
    # Print summary to console
    print("=== Distribution Analysis ===")
    print(f"Sample size: {stats_dict['n']}")
    print(f"Positive values: {stats_dict['n_positive']}")
    print(f"Mean: {stats_dict['mean']:.4f}")
    print(f"Std: {stats_dict['std']:.4f}")
    print(f"Coefficient of Variation (CV): {stats_dict['cv']:.4f}")
    print(f"CV > 1 (heavy-tailed indicator): {stats_dict['cv_heavy_tailed']}")
    
    if 'qq_correlation' in stats_dict and stats_dict['qq_correlation'] is not None:
        print(f"QQ-plot correlation with exponential: {stats_dict['qq_correlation']:.4f}")
        if stats_dict['qq_correlation'] < 0.95:
            print("  → Data likely NOT exponential")
    
    return fig, axes, stats_dict



def main():

    output_file = "./user_projects/output_sensibility_analysis/membrane_integrity_groups/results_summary.txt"
    X, Y = extract_X_Y(output_file, 'Group_m1')
    #plot_histogram_with_analysis(Y)
    plot_scatter_sets(Y, X, set_names=['Group_m1'])

    result_file = "./output_integrity3_groups/membrane_integrity.txt"
    param_names_file = "./output_integrity3_groups/param_names.txt"
    param_values_file = "./output_integrity3_groups/param_values_membrane_integrity.txt"
    
    #names = ["_".join(p) for p in param]
    groups = ["Group_epib1", "Group_epib2", "Group_epib2", "Group_m2", "Group_m2", "Group_m1", "Group_m3", "Group_m3"]

    #combine_files_with_header(result_file, param_values_file, output_file, groups)
    
    bounds = [[0.0, 1.0] * len(groups)]

    #Si = analyze_sobol(result_file, param_names_file, bounds, groups=groups)
    #Si.plot()
    #plt.show()
    
if __name__ == '__main__':
    main()