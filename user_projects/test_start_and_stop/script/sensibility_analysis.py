import os
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
            for cell_def in root.find(path_to_cell_def).iter("cell_definition"):
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
def membrane_integrity(cells_info, membrane_type, membrane_thickness, X_min=10000, X_max=10000):
    positions, cells_type = cells_info
    
    mask = cells_type == membrane_type
    membrane_positions = positions[mask]
    x_filter_membrane = (membrane_positions[:,0] >= X_min) & (membrane_positions[:,0] <= X_max)

    y_range_membrane = x_filter_membrane[:,1].max() - x_filter_membrane[:,1].min()

    if y_range_membrane > membrane_thickness:
        return False
    else:
        return True


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


def get_output_files(path, prefix='output', suffix='.mat'):
    
    pattern = os.path.join(path, prefix + '*' + suffix)
    files = glob.glob(pattern)

    files.sort()

    return files

def concatenate_results(files_mat, files_xml):
    return


def define_set_param(num_vars, names, bounds): 
    print("Defining parameter space and sampling...")
    problem = {
        'num_vars': num_vars,
        'names': names,
        'bounds': bounds
    }

    param_values = sample(problem, 2) #Génère N*(2+D) jeux de paramètres avec D le nombre de paramètres et N un multiple de 2 fourni en argument
    return problem, param_values

def get_physicell_output(param_values): 
    output_storage = np.zeros(param_values.shape[0]) 
    for i in range (len(param_values)): 
        print(f"Running simulation {i+1} with parameters: life_cycle_cancer_cell={param_values[i][0]}, transformation_rate_epi_basal_cancer={param_values[i][1]}")
        life_cycle_cancer_cell = param_values[i][0]
        transformation_rate_epi_basal_cancer = param_values[i][1]
        dict_xml = define_settings(life_cycle_cancer_cell, transformation_rate_epi_basal_cancer)

        xml_element = dict_to_xml(dict_xml)
        tree = ET.ElementTree(xml_element)
        tree.write(os.path.join(root_path, "PhysiCell/config/PhysiCell_settings.xml"), encoding="utf-8", xml_declaration=True) #path à changer

        process0 = subprocess.run(
        ["rm", "-rf", "*"],
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell/output") #to clean the output folder before each simulation
        )

        print("Configuration file updated. Starting simulation...")
        # Run .exe file 
        process1 = subprocess.run(
        [os.path.join(root_path, "PhysiCell/project.exe")], # !! To change if linux 
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell")
        )

        # Print the output and error (if any)
        print("Output:")
        print(process1.stdout)
        print("Error:")
        print(process1.stderr)

        #get output from output folder 
        output_path = os.path.join(root_path, "PhysiCell/output")  
        files_by_timestep = list_path_folder(output_path)
        result_mat = get_matrix_ids(files_by_timestep, root_path)
        area_over_time = computation_area_over_time(result_mat)
        output_storage[i] = area_over_time

        #stocker la vidéo 
        process2 = subprocess.run(
        ["make", "jpeg"],
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell")
        )
        process3 = subprocess.run(
        ["make", "gif"],
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell")
        )
        process4 = subprocess.run(
        ["make", "movie"],
        capture_output=True,
        text=True,
        cwd=os.path.join(root_path, "PhysiCell")
        )

        print("Output:")
        print(f"{process2.stdout}, {process3.stdout}, {process4.stdout}")
        print("Error:")
        print(f"{process2.stderr}, {process3.stderr}, {process4.stderr}")

        src = os.path.join(root_path, "PhysiCell/output/out.mp4")
        dst_folder = os.path.join(root_path, "output_video")
        dst = os.path.join(dst_folder, f"movie_simulation_{i+1}.mp4")
        os.makedirs(dst_folder, exist_ok=True)
        # Move video file and rename it
        shutil.move(src, dst)

    return output_storage

def main():
    #write the path of all the nodes we want to modify 
    #make a copy of the original .xml in the output folder
    #for a node in the xml
        #store original value

        #for all the value to test
            #modify the .xml functions works
            #run the simulation
            #perform analysis
            #store analysis results
            #erase useless output    
        #restore the value to original    
        
    #import argparse
    
    #parser = argparse.ArgumentParser(description='Convert PhysiCell output to ParaView format')
    #parser.add_argument('output_dir', help='Directory containing .mat and .xml output files')
    #parser.add_argument('--clean', action='store_true', help='Remove existing output files before processing')
    #parser.add_argument('--prefix', default='timestep', help='File prefix for VTU files (default: timestep)')
    
    #args = parser.parse_args()
    
    #pvd_file, vtu_files = create_pvd(args.output_dir, args.clean, args.prefix)
    
    #if pvd_file:
    #    print(f"\nConversion complete! To visualize these files in ParaView:")
    #    print(f"1. Open ParaView")
    #    print(f"2. File > Open > Navigate to: {os.path.abspath(pvd_file)}")
    #    print(f"3. Click 'Apply' in the Properties panel to load the data")
    print(os.getcwd())
    
    mat_file = "./output/output00000001_cells.mat"
    initial_xml_file = "./output/initial.xml"
    #print(parse_physicell_labels(xml_file))
    #print(analyse_output(mat_file, xml_file))

    path = "./cell_definitions/cell_definition/phenotype/cell_transformations/transformation_rates/transformation_rate"
    xml_file = "./config/PhysiCell_settings.xml"
    value = 0.0
    #print(parse_physicell_labels(initial_xml_file))
    #cells_info = extract_position_type_data(mat_file, initial_xml_file)
    #heights = cell_type_height(cells_info, -200, 200, 20, 2)
    #plot_cells_2D(heights)
    #print(cell_above([-50,-170], heights))
    #print(modify_xml(xml_file, path, value, name_cell_def="epi_basal", name_interact_cell_def="epi_inter1"))

    path = "./output/"
    #print(get_output_files(path, prefix='output', suffix='.mat'))
    print(define_set_param())

if __name__ == '__main__':
    main()