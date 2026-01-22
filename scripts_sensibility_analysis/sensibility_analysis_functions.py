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


def modify_xml(xml_file, path, value, name_cell_def="", name_interact_cell_def="", variable_name=""):
    """Modify XML file with error handling for incorrect paths and names"""
    try:
        if not os.path.exists(xml_file):
            raise FileNotFoundError(f"XML file '{xml_file}' does not exist")
        
        target_interaction = None
        subroot = None
        path_to_attribute = ""
        element = None

        modified = False
        tree = ET.parse(xml_file)
        root = tree.getroot()

        if variable_name != "":
            path_to_variable = path.split("variable/", 2)[0]
            path_to_attribute = "./" + path.split("variable/", 2)[1]
            
            if root.find(path_to_variable) is None:
                raise ValueError(f"Path: '{path_to_variable}' not found in XML structure")
            
            for var in root.findall(os.path.join(path_to_variable,"variable")):
                if var.get("name") == variable_name:
                    subroot = var
                    break

            if subroot is None:
                raise ValueError(f"Variable: '{variable_name}' not found in XML structure")

        if name_cell_def != "": 
            path_to_cell_def = path.split("cell_definition/", 2)[0]
            path_to_attribute = "./" + path.split("cell_definition/", 2)[1]
            
            if root.find(path_to_cell_def) is None:
                raise ValueError(f"Path: '{path_to_cell_def}' not found in XML structure")
            
            for cell_def in root.findall(os.path.join(path_to_cell_def,"cell_definition")):
                if cell_def.get("name") == name_cell_def:
                    subroot = cell_def
                    break
            if subroot is None:
                raise ValueError(f"Cell definition: '{name_cell_def}' not found in XML structure")
            
            if name_interact_cell_def != "":
                for interact_cell_def in subroot.findall(path_to_attribute):
                    if interact_cell_def.get("name") == name_interact_cell_def:
                        target_interaction = interact_cell_def
                        break
                
                if target_interaction is None:
                    error_msg = f"Interaction cell definition '{name_interact_cell_def}' not found in cell definition '{name_cell_def}' at path '{path_to_attribute}'"
                    raise ValueError(error_msg)
    
        if target_interaction:
            target_interaction.text = str(value)
            modified = True
        elif subroot != None and (path_to_attribute != ""):
            element = subroot.find(path_to_attribute)
            if element != None:
                element.text = str(value)
                modified = True
            else:
                raise ValueError("Invalid path in the xml")


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

def analyze_sobol(result_file, param_names_file, bounds, groups=[], column = 1):
    import ast
    with open(param_names_file, 'r') as f:
        names = [ast.literal_eval(line.strip()) for line in f if line.strip()]
    names = np.asarray([f"{';'.join(sublist)}" for sublist in names])

    problem = define_problem(len(names), names, bounds, groups)
    print(problem)
    
    with open(result_file, 'r') as f:
        lines = f.readlines()

    sorted_lines = sorted(lines, key=lambda x: int(x.split()[0]))
    Y = np.asarray([float(line.split()[column]) for line in sorted_lines])
    print("Array:", Y)
    print("Length:", len(Y))
    
    Si = analyze(problem, Y, print_to_console=True)
    return Si



###########  Functions for output processing and analysis ###########

def process_file(input_file, output_file):
    """
    Process the file to remove duplicate entries where first column is identical
    AND second column values are the same.
    """
    data_dict = {}

    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    for line_num, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
            
        parts = line.split()
        if len(parts) >= 2:
            first_col = parts[0]
            second_col = parts[1]
            
            if first_col not in data_dict:
                data_dict[first_col] = []
            data_dict[first_col].append((second_col, line))
    
    output_lines = []
    
    for first_col in sorted(data_dict.keys()):
        entries = data_dict[first_col]
        
        if len(entries) == 1:
            output_lines.append(entries[0][1])
        else:
            first_value = entries[0][0]
            all_same = all(entry[0] == first_value for entry in entries[1:])
            
            if all_same:
                output_lines.append(entries[0][1])
                print(f"Removed {len(entries)-1} duplicate(s) for {first_col} (all values = {first_value})")
            else:
                for entry in entries:
                    output_lines.append(entry[1])
                print(f"WARNING: Different values for {first_col}: {[e[0] for e in entries]}")
    
    with open(output_file, 'w') as f:
        for line in output_lines:
            f.write(line + '\n')
    
    original_count = sum(1 for line in lines if line.strip())
    new_count = len(output_lines)
    
    print(f"\nResults:")
    print(f"Original file: {original_count} lines")
    print(f"Processed file: {new_count} lines")
    print(f"Lines removed: {original_count - new_count}")
    print(f"Output saved to: {output_file}")
    
    return new_count


def combine_files_with_header(file1_path, file2_path, output_path, column_names, column_1='nbr_breaks'):
    
    df1 = pd.read_csv(file1_path, sep=r'\s+', header=None, names=['row_idx', column_1])
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
            nbr_breaks_list[row_idx] = row[column_1]
        
    df_combined = pd.DataFrame({column_1: nbr_breaks_list})
    df_combined = pd.concat([df_combined, df2.reset_index(drop=True)], axis=1)
        
    df_combined.to_csv(output_path, sep='\t', index=False, float_format='%.6e')
        
    print(f"Combined file created: {output_path}")
    print(f"Shape of combined data: {df_combined.shape}")
    print(f"Rows from file1 assigned: {len(df1[df1['row_idx'] < len(df2)])}")
        
    return df_combined

def combine_files_with_header_2(file1_path, file2_path, output_path, column_names, column_1, column_2):
    
    df1 = pd.read_csv(file1_path, sep=r'\s+', header=None, names=['row_idx', column_1, column_2])
    df2 = pd.read_csv(file2_path, sep=r'\s+', header=None)

    print(df1)
        
    if len(column_names) == df2.shape[1]:
        df2.columns = column_names
    
    else:
        print(f"Warning: column_names has {len(column_names)} items, but file2 has {df2.shape[1]} columns")
        df2.columns = [f"col_{i+1}" for i in range(df2.shape[1])]
        
    value_column1_list = [0.0] * len(df2)
    value_common2_list = [0.0] * len(df2)
        
    for idx, row in df1.iterrows():
        row_idx = int(row['row_idx'])  
        if 0 <= row_idx < len(df2):  
            value_column1_list[row_idx] = row[column_1]
            value_common2_list[row_idx] = row[column_2]

        
    df_combined = pd.DataFrame({column_1: value_column1_list, column_2:value_common2_list })
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

def save_dataframes_to_txt(df_list, filename, 
                          df_names=None, include_index=True):
    with open(filename, 'w') as f:
        for i, df in enumerate(df_list):
            if df_names and i < len(df_names):
                f.write(f"=== DataFrame: {df_names[i]} ===\n")
            else:
                f.write(f"=== DataFrame {i+1} ===\n")
            
            df_str = df.to_string(index=include_index, header=True)
            f.write(df_str)
            f.write("\n\n")  
    
    print(f"Saved {len(df_list)} DataFrames to '{filename}'")


def load_dataframes_from_txt(filename):
    dataframes = []
    current_df_lines = []
    reading_df = False
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        # Check if DataFrame header
        if line.startswith("=== DataFrame"):
            if reading_df and current_df_lines:
                df_str = ''.join(current_df_lines)
                try:
                    lines_split = df_str.strip().split('\n')
                    if len(lines_split) > 1:
                        df = pd.read_csv(pd.io.common.StringIO(df_str), 
                                        sep=r'\s+', engine='python')
                        dataframes.append(df)
                except:
                    print(f"Warning: Could not parse DataFrame section")
                
                current_df_lines = []
            
            reading_df = True
            continue
        
        if reading_df and line.strip():
            current_df_lines.append(line)
    
    if reading_df and current_df_lines:
        df_str = ''.join(current_df_lines)
        try:
            df = pd.read_csv(pd.io.common.StringIO(df_str), 
                            sep=r'\s+', engine='python')
            dataframes.append(df)
        except:
            print(f"Warning: Could not parse last DataFrame section")
    
    return dataframes

