from sensibility_analysis_functions import *
from saved_data_modifications import *
import subprocess


# Write a function to compute the thickness of the epi inter, and its growth rate
def compute_epi_thickness(mat_file, xml_file):
    positions, types, cells_ID = extract_position_type_data(mat_file, xml_file)

    mask_epi_inter_cells = types == 1.0
    positions_epi_inter = positions[mask_epi_inter_cells]
    positions_y_epi_inter = np.sort(positions_epi_inter[:, 1])

    min = np.average(positions_y_epi_inter[:10])
    max = np.average(positions_y_epi_inter[-10:])

    return abs(max - min)


# Write a function to test if the epithelium is stable in size
def compute_epi_stability(mat_files, xml_file, sort=False):
    if sort:
        # Sorting only works if your files are named acccording to PhysiCell standard outputs
        mat_files.sort(key=lambda x: int(x.split("_")[0].split("output")[-1]))
        print(mat_files)

    sizes = []
    for file in mat_files:
        sizes.append(compute_epi_thickness(file, xml_file))

    growth_rates = []
    for i in range(1, len(sizes)):
        growth_rates.append(sizes[i] - sizes[i - 1])
    return np.average(growth_rates)


# Write a function to track the number of the cells type
def compute_number_cells(mat_file, xml_file, tracked_types):
    positions, types, cells_ID = extract_position_type_data(mat_file, xml_file)
    types_dict = {}

    for t in tracked_types:
        mask = types == t
        types_dict[t] = [len(types[mask])]

    return types_dict


def compute_number_cells_over_time(mat_files, xml_file, tracked_types):
    pop_dict = {}
    for t in tracked_types:
        pop_dict[t] = []

    for file in mat_files:
        dict = compute_number_cells(file, xml_file, tracked_types)
        for type in tracked_types:
            pop_dict[type] = pop_dict[type] + dict[type]

    return pop_dict


# Write a function to evaluate the membrane thickness


def modify_csv(file_path, cell_rule, value):
    # Validate cell_rule format
    if len(cell_rule) != 3 or cell_rule[0] != "cell_rule":
        raise ValueError(
            "cell_rule must be in format: ['cell_rule', row_index, col_index]"
        )

    # Get row and column indices
    row_idx = cell_rule[1]
    col_idx = cell_rule[2]

    # Read all lines from the file
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Check if row index is valid
    if row_idx < 0 or row_idx >= len(lines):
        raise IndexError(
            f"Row index {cell_rule[1]} is out of range. File has {len(lines)} rows."
        )

    # Split the line by commas and modify the specific column
    line_parts = lines[row_idx].strip().split(",")

    # Check if column index is valid
    if col_idx < 0 or col_idx >= len(line_parts):
        raise IndexError(
            f"Column index {cell_rule[2]} is out of range. Row {cell_rule[1]} has {len(line_parts)} columns."
        )

    # Modify the specific cell
    line_parts[col_idx] = str(value)

    # Reconstruct the line
    lines[row_idx] = ",".join(line_parts) + "\n"

    # Create the output file path
    folder = os.path.dirname(file_path)
    output_path = os.path.join(folder, "temp_cell_rules.csv")

    # Write all lines to the new file
    with open(output_path, "w") as file:
        file.writelines(lines)

    print(f"Modified copy saved to: {output_path}")
    print(
        f"Modified cell at row {cell_rule[1]}, column {cell_rule[2]} to value: {value}"
    )

    return output_path


def evaluate_epi_growth(
    xml_file,
    param_treepaths,
    param_values,
    save_output_folder,
    temp_output,
    restart=False,
    restart_int=0,
    end_int=0,
):
    output_folder = os.path.join(os.getcwd(), save_output_folder)
    save_output = os.path.join(output_folder, "run_output.txt")
    output_storage_file = os.path.join(output_folder, "epi_growth.txt")
    output_storage_cell_pop = os.path.join(output_folder, "epi_cell_pop.txt")
    temp_output_folder = os.path.join(os.getcwd(), temp_output)
    temp_xml_file = os.path.join(temp_output_folder, os.path.basename(xml_file))
    start_file = 0

    # Identify the cell_rules file if defined
    tree = ET.parse(xml_file)
    root = tree.getroot()
    filepath = ""
    if root.find("./cell_rules/rulesets/ruleset/folder") != None:
        folder = root.find("./cell_rules/rulesets/ruleset/folder").text
        file = root.find("./cell_rules/rulesets/ruleset/filename").text
        print(file)
        print(folder)
        filepath = os.path.join(folder, file)
        print(filepath)

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
            with open(os.path.join(output_folder, "param_names.txt"), "w") as fp:
                for item in param_treepaths:
                    fp.write("%s\n" % item)
            np.savetxt(os.path.join(output_folder, "param_values.txt"), param_values)

        if os.path.exists(xml_file):
            shutil.copy(xml_file, temp_xml_file)

        start_file = restart_int
    else:
        if os.path.exists(xml_file) & (not os.path.exists(temp_xml_file)):
            shutil.copy(xml_file, temp_xml_file)

        param_values = np.loadtxt(os.path.join(output_folder, "param_values.txt"))
        pattern = r'["\'](.*?)["\']'
        with open(
            os.path.join(os.path.join(output_folder, "param_names.txt")), "r"
        ) as file:
            param_treepaths = [
                re.findall(pattern, line.strip()) for line in file if line.strip()
            ]
        if restart_int == 0:
            # get_output_files returns a list of files with corresponding prefix and suffix in folder sorted by the run number in their name
            files = get_output_files(output_folder, prefix="out_", suffix=".gif")
            last_file = os.path.basename(files[-1])
            last_file_int = int(last_file[4 : len(last_file) - 4])
            start_file = last_file_int + 1
        else:
            start_file = restart_int
    if end_file < start_file:
        print(
            f"The end index({end_file}) is inferior to the start index,({start_file}) check your input values"
        )
        raise ValueError

    if not os.path.isfile(save_output):
        with open(save_output, "w") as f:
            pass
    if not os.path.isfile(output_storage_file):
        with open(output_storage_file, "w") as f:
            pass
    if not os.path.isfile(output_storage_cell_pop):
        with open(output_storage_cell_pop, "w") as f:
            pass

    output_growth_rates = []
    output_epi_sizes = []
    output_cell_pop = []
    # Changing output folder
    modify_xml(temp_xml_file, "save/folder", temp_output)
    '''
    # Create a copy and modify the start and stop data folder
    shutil.copytree(
        "./config/start_and_stop_saving_files",
        os.path.join(temp_output_folder, "start_and_stop_saving_files",),
        dirs_exist_ok=True
    )
    modify_xml(
        temp_xml_file,
        "user_parameters/saving_folder",
        os.path.join(temp_output_folder, "start_and_stop_saving_files"),
    )
    modify_xml(
        temp_xml_file,
        "user_parameters/init_cells_filename",
        os.path.join(temp_output_folder, "start_and_stop_saving_files/initial.tsv"),
    )
    '''

    for i in range(end_file - start_file):
        print(" param values: ",param_values[i + start_file, :])
        print("param treepaths: ",param_treepaths)
        for j, val in enumerate(param_values[i + start_file, :]):
            print("j = " ,j, " val = ", val)
            name_cell_def = ""
            name_interact_cell_def = ""
            variable_name = ""
            substrate = ""

            if len(param_treepaths[j]) > 1:
                name_cell_def = param_treepaths[j][1]
            if len(param_treepaths[j]) > 2:
                name_interact_cell_def = param_treepaths[j][2]
            if len(param_treepaths[j]) > 3:
                variable_name = param_treepaths[j][3]
            if len(param_treepaths[j]) > 4:
                substrate = param_treepaths[j][4]

            if param_treepaths[j][0] == "cell_rule":
                modify_csv(filepath, param_treepaths[j], val)
                modify_xml(
                    temp_xml_file,
                    "./cell_rules/rulesets/ruleset/filename",
                    "temp_cell_rules.csv",
                )
            else:
                modify_xml(
                    temp_xml_file,
                    param_treepaths[j][0],
                    val,
                    name_cell_def=name_cell_def,
                    name_interact_cell_def=name_interact_cell_def,
                    variable_name=variable_name,
                    substrate=substrate,
                )
            # After modifying the .xml we need to modify the saved start_and_stop files for proper initialization
            # modify_cell_data(os.path.join(temp_output_folder, 'start_and_stop_saving_files'),)
        # Running simulation
        process0 = subprocess.run(
            ["./test_death", temp_xml_file], capture_output=True, text=True
        )

        # Outputs analysis
        label_file = os.path.join(temp_output_folder, "initial.xml")
        mat_files = get_output_files(temp_output_folder)
        # Average growth rate of epi_inter layer over the run
        output_growth_rates.append(
            np.average(compute_epi_stability(mat_files, label_file))
        )
        # Size of the epi_inter layer at the last time step
        output_epi_sizes.append(compute_epi_thickness(mat_files[-1], label_file))
        # Cells population size
        output_cell_pop.append(
            compute_number_cells_over_time(mat_files, label_file, [0, 1])
        )
        print("Run ", i + start_file, " completed")

        process1 = subprocess.run(
            ["make", "gif", f"OUTPUT={temp_output_folder}"],
            capture_output=True,
            text=True,
        )

        shutil.copyfile(
            f"{temp_output_folder}/out.gif", f"{output_folder}/out_{i + start_file}.gif"
        )
        with open(output_storage_file, "a") as f:
            f.write(
                f"{i + start_file} {output_growth_rates[i]} {output_epi_sizes[i]}\n"
            )
        with open(output_storage_cell_pop, "a") as f:
            f.write(f"{i + start_file} {output_cell_pop[i]}\n")
        with open(save_output, "a") as f:
            f.write(process0.stdout)
        with open(save_output, "a") as f:
            f.write(process1.stdout)

    # delete temp output at the end of the run
    if os.path.isdir(temp_output_folder):
        shutil.rmtree(temp_output_folder)

    return output_growth_rates, output_epi_sizes


def main():
    ### Analysis for membrane_growth_run_2.py results ###
    """
    param = [["cell_definitions/cell_definition/phenotype/secretion/substrate/secretion_rate", "epi_inter", "", "", "div_inhib"],
         ["microenvironment_setup/variable/physical_parameter_set/diffusion_coefficient" , "", "", "div_inhib"],
         ["microenvironment_setup/variable/physical_parameter_set/decay_rate", "", "", "div_inhib"]]

    param_bounds = [[1.0, 1000.0],
          [500.0, 1500.0],
          [0.05, 0.5]]


    names = ["_".join(p) for p in param]
    num_vars = len(param)
    bounds = param_bounds
    xml_file = "./config/PhysiCell_settings.xml"
    param_values = define_set_param(num_vars, names, bounds, sample_size=128)

    total_param = param
    
    #Saving sensibility analysis sobol results groups
    result_file = "./output_sensibility_analysis/output_growth_2/epi_growth.txt"
    #process_file(result_file, result_file)
    param_names_file = "./output_growth_2/output_growth/param_names.txt"
    param_values_file = "./output_growth_2/output_growth/param_values.txt"

    #Si = analyze_sobol(result_file, param_names_file, bounds, column=2)
    #save_dataframes_to_txt([Si.to_df()[0], Si.to_df()[1], Si.to_df()[2]], "./output_sensibility_analysis/output_growth_2/data_output_growth_2.txt")

    #Si.plot()
    #plt.show()

    #Combine file with header
    output_path = "./output_sensibility_analysis/output_growth_2/result_summary_growth_2.txt"
    #combine_files_with_header_2(result_file, param_values_file, output_path, names, 'growth_rate', 'epi_size')
    
    #Plotting scatter plot groupped variable 
    X, Y = extract_X_Y(output_path, names[2], Y_column='growth_rate')
    plot_scatter_sets(Y, X, set_names="", xlabel="secretion rate", title ="Growth rate in function of secretion rate")
    """

    ### Analysis for membrane_growth_run.py results ###
    param = [["cell_definitions/cell_definition/phenotype/secretion/substrate/secretion_rate", "epi_inter", "", "", "div_inhib"],
         ["microenvironment_setup/variable/physical_parameter_set/diffusion_coefficient" , "", "", "div_inhib"],
         ["microenvironment_setup/variable/physical_parameter_set/decay_rate", "", "", "div_inhib"]]

    param_bounds = [[1.0, 1000.0],
          [500.0, 1500.0],
          [0.05, 0.5]]

    cell_rules = [['cell_rule',1, 5]]
    cell_rules_bounds = [[0.01, 4]]

    names = ["_".join(p) for p in param] + ["_".join(str(r)) for r in cell_rules]
    num_vars = len(param) + len(cell_rules)
    bounds = param_bounds + cell_rules_bounds
    xml_file = "./config/PhysiCell_settings.xml"
    param_values = define_set_param(num_vars, names, bounds, sample_size=128)
    print("Param values shape: ", param_values.shape)

    total_param = param + cell_rules    

    ### Saving sensibility analysis sobol results groups ###
    result_file = "./output_sensibility_analysis/output_growth/epi_growth.txt"
    process_file(result_file, result_file)
    param_names_file = "./output_growth/param_names.txt"
    param_values_file = "./output_growth/param_values.txt"

    #Combine file with header
    output_path = "./output_sensibility_analysis/output_growth/result_summary_growth.txt"
    #combine_files_with_header_2(result_file, param_values_file, output_path, names, 'growth_rate', 'epi_size')


    #Si = analyze_sobol(result_file, param_names_file, bounds, column=2)
    #save_dataframes_to_txt([Si.to_df()[0], Si.to_df()[1], Si.to_df()[2]], "./output_sensibility_analysis/output_growth/data_output_growth.txt")
    #Si.plot()
    #plt.show()

    #Plotting scatter plot groupped variable 
    X, Y = extract_X_Y(output_path, names[3], Y_column='growth_rate')
    plot_scatter_sets(Y, X, set_names="", xlabel="Division downregulation", title ="Epithelium growth rate in function of division downregulation")

if __name__ == "__main__":
    main()
