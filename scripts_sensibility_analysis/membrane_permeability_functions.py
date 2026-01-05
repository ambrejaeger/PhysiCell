from sensibility_analysis_functions import *

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
        print("This runs")
        if os.path.exists(xml_file) & (not os.path.exists(temp_xml_file)):
            shutil.copy(xml_file, temp_xml_file)

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

def main():
    param = [["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "attracted"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "attracted", "membrane"],
    ["cell_definitions/cell_definition/phenotype/motility/migration_bias", "attracted"],
    ["cell_definitions/cell_definition/phenotype/volume/total", "attracted"]
    ]

    bounds = [[0.0, 1.0],
          [0.0, 1.0],
          [0.0, 1.0],
          [0.0, 1.0],
          [3000, 5500]]
    """
    names = ["_".join(p) for p in param]
    groups = ['Group_m1', 'Group_m2', 'Group_m2', 'Group_attr', 'Group_vattr']
    result_file = "./output_integrity4_groups/membrane_integrity.txt"
    param_names_file = "./output_integrity4_groups/param_names.txt"
    param_values_file = "./output_integrity4_groups/param_values_membrane_permeability.txt"
    output_file = "./output_sensibility_analysis/membrane_permeability_groups/result_file.txt"
    Si = analyze_sobol(result_file, param_names_file, bounds, groups=groups)
    save_dataframes_to_txt([Si.to_df()[0], Si.to_df()[1], Si.to_df()[2]], "./output_sensibility_analysis/membrane_permeability_groups/data_sensibility_analysis.txt")
    """

    # OR FOR NO GROUPS ANALYSIS
    """
    result_file = "./output_integrity4/membrane_integrity.txt"
    param_names_file = "./output_integrity4/param_names.txt"
    Si = analyze_sobol(result_file, param_names_file, bounds)
    save_dataframes_to_txt([Si.to_df()[0], Si.to_df()[1], Si.to_df()[2]], "./output_sensibility_analysis/membrane_permeability_no_group/data_sensibility_analysis.txt")
    """
    
    #Si.plot()
    #plt.show()

    ###### Obtaining Histogram plot permeability groups ######
    '''combine_files_with_header(result_file, param_values_file, output_file, groups, column_1="cross")
    output_file = "./output_sensibility_analysis/membrane_permeability_groups/result_file.txt" 
    X, Y = extract_X_Y(output_file, 'Group_m1', Y_column="cross")
    plt.hist(X, bins=20, weights = Y)
    plt.title("Succesful crossing in funtion of adhesion affinity")
    plt.xlabel("membrane membrane adhesion affinity")
    plt.ylabel("Runs")
    plt.show()'''

    ###### Obtaining Histogram plot permeability NO groups ######
    """result_file = "./output_integrity4/membrane_integrity.txt"
    param_values_file = "./output_integrity4/param_values_membrane_permeability.txt"
    output_file = "./output_sensibility_analysis/output_permeability_no_group/result_file.txt" 
    combine_files_with_header(result_file, param_values_file, output_file, names, column_1="cross")
    
    X, Y = extract_X_Y(output_file, 'cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity_membrane_membrane', Y_column="cross")
    plt.hist(X, bins=20, weights = Y)
    plt.title("Succesful crossing in funtion of adhesion affinity")
    plt.xlabel("membrane membrane adhesion affinity")
    plt.ylabel("Runs")
    plt.show()"""

    #plot_scatter_sets(Y, X, set_names=['Group_m1'], xlabel="adhesion affinity membrane membrane", title="Cross in function of adhesion affinity")
    
if __name__ == "__main__":
    main()
