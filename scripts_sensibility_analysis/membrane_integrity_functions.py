from sensibility_analysis_functions import *


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


def compute_membrane_integrity(mat_files, temp_output_folder, membrane_type, tolerance, X_min=-10000, X_max=10000):
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
        output_storage[i] = math.fsum(compute_membrane_integrity(mat_files, temp_output_folder, 2, tolerance))
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
    groups = ["Group_epib1", "Group_epib2", "Group_epib2", "Group_m2", "Group_m2", "Group_m1", "Group_m3", "Group_m3"]
    
    """
    #Saving sensibility analysis sobol results no group
    result_file = "./output_integrity3/membrane_integrity.txt"
    param_names_file = "./output_integrity3/param_names.txt"
    Si = analyze_sobol(result_file, param_names_file, bounds)
    save_dataframes_to_txt([Si.to_df()[0], Si.to_df()[1], Si.to_df()[2]], "./output_sensibility_analysis/membrane_integrity_no_group/data_sensibility_analysis.txt")
    """

    #Saving sensibility analysis sobol results groups
    result_file = "./output_integrity3_groups/membrane_integrity.txt"
    param_names_file = "./output_integrity3_groups/param_names.txt"
    Si = analyze_sobol(result_file, param_names_file, bounds, groups = groups)
    save_dataframes_to_txt([Si.to_df()[0], Si.to_df()[1], Si.to_df()[2]], "./output_sensibility_analysis/membrane_integrity_groups/data_sensibility_analysis.txt")
    

    output_file = "./output_integrity3/membrane_integrity_processed.txt"    
    result_file = "./output_integrity3/membrane_integrity.txt"

    result_file = "./output_integrity3_groups/membrane_integrity.txt"
    param_names_file = "./output_integrity3_groups/param_names.txt"
    param_values_file = "./output_integrity3_groups/param_values_membrane_integrity.txt"
    
   

    #Si = analyze_sobol(result_file, param_names_file, bounds)

    #OR FOR GROUPPED VARIABLES
    #Si = analyze_sobol(result_file, param_names_file, bounds, groups=groups)

    #Si.plot()
    #plt.show()

    #output_file = "./user_projects/output_sensibility_analysis/membrane_integrity_groups/results_summary.txt"
    #process_file(result_file, output_file )

    """
    #Plotting scatter plot groupped variable
    output_file = "./output_sensibility_analysis/membrane_integrity_groups/results_summary.txt" 
    X, Y = extract_X_Y(output_file, 'Group_m1')
    plot_scatter_sets(Y, X, set_names=['Group_m1'], xlabel="adhesion affinity membrane membrane", title="Number of breaks in function of adhesion affinity")
    """

if __name__ == '__main__':
    main()