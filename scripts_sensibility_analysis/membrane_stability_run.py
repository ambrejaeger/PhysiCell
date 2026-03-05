from membrane_growth_functions import *
from membrane_apoptosis_functions import *
from saved_data_modifications import *
import argparse

param = [
    [
        "cell_definitions/cell_definition/phenotype/secretion/substrate/secretion_rate",
        "epi_inter",
        "",
        "",
        "div_inhib",
    ],
    [
        "microenvironment_setup/variable/physical_parameter_set/diffusion_coefficient",
        "",
        "",
        "div_inhib",
    ],
    [
        "microenvironment_setup/variable/physical_parameter_set/decay_rate",
        "",
        "",
        "div_inhib",
    ],
    [
        "cell_definitions/cell_definition/phenotype/death/model/death_rate",
        "epi_inter",
        "",
        "",
        "",
    ]
]

saved_data_path = [["Secretion/Secretion 1/Secretion_Rate", "type", 1], [], [], ["Death/Model 0/Rate", "type", 1]]

param_bounds = [[1.0, 1000.0], [500.0, 1500.0], [0.05, 100], [0.000001, 0.001]]

cell_rules = [['cell_rule',"1", "6"], ['cell_rule',"2", "6"]]
cell_rules_bounds = [[0.01, 4], [0.01, 4]]

names = ["_".join(p) for p in param] + ["_".join(str(r)) for r in cell_rules]
num_vars = len(param) + len(cell_rules)
bounds = param_bounds + cell_rules_bounds
xml_file = "./config/PhysiCell_settings.xml"
param_values = define_set_param(num_vars, names, bounds, sample_size=64)

total_param = param + cell_rules

''' 
def evaluate_epi_growth_2(
    xml_file,
    param_treepaths,
    param_values,
    save_output_folder,
    temp_output,
    restart=False,
    restart_int=0,
    end_int=0,
    saved_data_paths=[],
):
    output_folder = os.path.join(os.getcwd(), save_output_folder)
    save_output = os.path.join(output_folder, "run_output.txt")
    output_storage_file = os.path.join(output_folder, "epi_growth.txt")
    output_storage_cell_pop = os.path.join(output_folder, "epi_cell_pop.txt")
    temp_output_folder = os.path.join(os.getcwd(), temp_output)
    temp_xml_file = os.path.join(temp_output_folder, os.path.basename(xml_file))
    start_file = 0
    print(param_treepaths)
    if not os.path.isdir(temp_output_folder):
        os.makedirs(temp_output_folder, exist_ok=False)

    # Identify the cell_rules file if defined
    tree = ET.parse(xml_file)
    root = tree.getroot()
    filepath = ""
    cell_rule_file = ""
    if root.find("./cell_rules/rulesets/ruleset/folder") != None:
        folder = root.find("./cell_rules/rulesets/ruleset/folder").text
        file = root.find("./cell_rules/rulesets/ruleset/filename").text
        filepath = os.path.join(folder, file)

        # Creating a copy in the temp_output_folder
        cell_rule_file = shutil.copy(
            filepath, os.path.join(temp_output_folder, "cell_rules.csv")
        )
        print(cell_rule_file)

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

    # Create a copy and modify the start and stop data folder
    print(
        shutil.copytree(
            "./config/start_and_stop_saving_files",
            os.path.join(temp_output_folder, "start_and_stop_saving_files"),
            dirs_exist_ok=True,
        )
    )
    modify_xml(
        temp_xml_file,
        "user_parameters/saving_folder",
        os.path.join(temp_output_folder, "start_and_stop_saving_files/"),
    )
    modify_xml(
        temp_xml_file,
        "user_parameters/init_cells_filename",
        os.path.join(temp_output_folder, "start_and_stop_saving_files/initial.tsv"),
    )

    for i in range(end_file - start_file):
        for j, val in enumerate(param_values[i + start_file, :]):
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
                print(param_treepaths)
                modify_csv(cell_rule_file, param_treepaths[j], val)
                modify_xml(
                    temp_xml_file,
                    "./cell_rules/rulesets/ruleset/folder",
                    temp_output_folder,
                )
                modify_xml(
                    temp_xml_file,
                    "./cell_rules/rulesets/ruleset/filename",
                    "cell_rules.csv",
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
            if param_treepaths[j][0].find("cell_definition") != -1:
                #cell_data_file = os.path.join(
                 #   temp_output_folder, "start_and_stop_saving_files/cell_data.txt"
                #)
                cell_data_file = os.path.join(temp_output_folder, "start_and_stop_saving_files/cell_data.txt")
                path = saved_data_paths[j]
                print(path[1])
                print(path[2])
                print(path[0])
                print(modify_cell_data(cell_data_file, path[1], path[2], path[0], val))
        
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
        
        dest_fin_svg = os.path.join(output_folder, f"final_{i + start_file}.svg")
        dest_init_svg = os.path.join(output_folder, f"initial_{i + start_file}.svg")
        initial_svg = os.path.join(temp_output_folder, "initial.svg")
        final_svg = os.path.join(temp_output_folder, "final.svg")
        shutil.copyfile(initial_svg, dest_init_svg)
        shutil.copyfile(final_svg, dest_fin_svg)

        with open(output_storage_file, "a") as f:
            f.write(
                f"{i + start_file} {output_growth_rates[i]} {output_epi_sizes[i]}\n"
            )
        with open(output_storage_cell_pop, "a") as f:
            f.write(f"{i + start_file} {output_cell_pop[i]}\n")
        with open(save_output, "a") as f:
            f.write(process0.stdout)
        #with open(save_output, "a") as f:
            #f.write(process1.stdout)

    # delete temp output at the end of the run
    if os.path.isdir(temp_output_folder):
        shutil.rmtree(temp_output_folder)

    return output_growth_rates, output_epi_sizes
'''

def evaluate_stability(
    xml_file,
    param_treepaths,
    param_values,
    save_output_folder,
    temp_output,
    restart=False,
    restart_int=0,
    end_int=0,
    saved_data_paths=[],
):
    output_folder = os.path.join(os.getcwd(), save_output_folder)
    save_output = os.path.join(output_folder, "run_output.txt")
    output_storage_file = os.path.join(output_folder, "epi_growth.txt")
    output_storage_cell_pop = os.path.join(output_folder, "epi_cell_pop.txt")
    temp_output_folder = os.path.join(os.getcwd(), temp_output)
    temp_xml_file = os.path.join(temp_output_folder, "temp_PhysiCell_settings.xml")
    start_file = 0
    print(param_treepaths)
    if not os.path.isdir(temp_output_folder):
        os.makedirs(temp_output_folder, exist_ok=False)

    #Make a copy of the file for the run
    # Creating a copy in the temp_output_folder
    shutil.copy(
        xml_file, temp_xml_file
    )

    # Identify the cell_rules file if defined
    tree = ET.parse(xml_file)
    root = tree.getroot()
    filepath = ""
    cell_rule_file = ""
    if root.find("./cell_rules/rulesets/ruleset/folder") != None:
        folder = root.find("./cell_rules/rulesets/ruleset/folder").text
        file = root.find("./cell_rules/rulesets/ruleset/filename").text
        filepath = os.path.join(folder, file)

        # Creating a copy in the temp_output_folder
        cell_rule_file = shutil.copy(
            filepath, os.path.join(temp_output_folder, "temp_cell_rules.csv")
        )
        print("cell rule folder modified: ",modify_xml(
                    temp_xml_file,
                    "./cell_rules/rulesets/ruleset/folder",
                    temp_output_folder,
        ))
        print("cell rule filename modified: ",modify_xml(
                    temp_xml_file,
                    "./cell_rules/rulesets/ruleset/filename",
                    "temp_cell_rules.csv",
        ))
        print(cell_rule_file)

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

        start_file = restart_int

    else:
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

    output_apoptosis_pos = []
    output_apoptosis_nbr = []
    output_growth_rates = []
    output_epi_sizes = []
    output_cell_pop = []
    # Changing output folder
    modify_xml(temp_xml_file, "save/folder", temp_output)

    # Create a copy and modify the start and stop data folder
    print(
        shutil.copytree(
            "./config/start_and_stop_saving_files",
            os.path.join(temp_output_folder, "start_and_stop_saving_files"),
            dirs_exist_ok=True,
        )
    )
    modify_xml(
        temp_xml_file,
        "user_parameters/saving_folder",
        os.path.join(temp_output_folder, "start_and_stop_saving_files/"),
    )
    modify_xml(
        temp_xml_file,
        "user_parameters/init_cells_filename",
        os.path.join(temp_output_folder, "start_and_stop_saving_files/initial.tsv"),
    )

    for i in range(end_file - start_file):
        for j, val in enumerate(param_values[i + start_file, :]):
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
                print(param_treepaths)
                modify_csv(cell_rule_file, param_treepaths[j], val)
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
            if param_treepaths[j][0].find("cell_definition") != -1:
                #cell_data_file = os.path.join(
                 #   temp_output_folder, "start_and_stop_saving_files/cell_data.txt"
                #)
                cell_data_file = os.path.join(temp_output_folder, "start_and_stop_saving_files/cell_data.txt")
                path = saved_data_paths[j]
                print(path[1])
                print(path[2])
                print(path[0])
                print(modify_cell_data(cell_data_file, path[1], path[2], path[0], val))
        
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
        # Average growth rate of epi_inter layer over the run
        output_apoptosis_pos.append(
            compute_median_apoptosis_position(label_file, mat_files)[1]
        )

        output_apoptosis_nbr.append(
            compute_median_apoptosis_position(label_file, mat_files)[0]
        )
        # Size of the epi_inter layer at the last time step
        output_epi_sizes.append(compute_epi_thickness(mat_files[-1], label_file))
        # Cells population size
        output_cell_pop.append(
            compute_number_cells_over_time(mat_files, label_file, [0, 1])
        )
        print("Run ", i + start_file, " completed")
        
        dest_fin_svg = os.path.join(output_folder, f"final_{i + start_file}.svg")
        dest_init_svg = os.path.join(output_folder, f"initial_{i + start_file}.svg")
        initial_svg = os.path.join(temp_output_folder, "initial.svg")
        final_svg = os.path.join(temp_output_folder, "final.svg")
        shutil.copyfile(initial_svg, dest_init_svg)
        shutil.copyfile(final_svg, dest_fin_svg)

        with open(output_storage_file, "a") as f:
            f.write(
                f"{i + start_file} {output_growth_rates[i]} {output_apoptosis_nbr[i]} {output_apoptosis_pos[i]} {output_epi_sizes[i]}\n"
            )
        with open(output_storage_cell_pop, "a") as f:
            f.write(f"{i + start_file} {output_cell_pop[i]}\n")
        with open(save_output, "a") as f:
            f.write(process0.stdout)
        #with open(save_output, "a") as f:
            #f.write(process1.stdout)

    # delete temp output at the end of the run
    if os.path.isdir(temp_output_folder):
        shutil.rmtree(temp_output_folder)

    return output_apoptosis_pos, output_epi_sizes

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script that runs ..._run.py file")
    parser.add_argument("--output", required=True, type=str)
    parser.add_argument("--temp_output", required=True, type=str)
    parser.add_argument(
        "--restart", required=False, type=str, choices=["0", "1", "true", "false"]
    )
    parser.add_argument("--start", required=False, type=int)
    parser.add_argument("--end", required=False, type=int)

    args = parser.parse_args()

    output = args.output
    temp_output = args.temp_output

    if args.restart and args.restart.lower() in ["1", "true"]:
        restart = True
    else:
        restart = False
    if args.start:
        start = args.start
    else:
        start = 0
    if args.end:
        end = args.end
    else:
        end = 0

    completed = False
    while not completed:
        if os.path.isdir(os.path.join(os.getcwd(), temp_output)):
            print(
                temp_output,
                " folder already exist in your working directory, by running this script you might overwrite an ongoing simulation.",
            )
            choice = input("Do you want to run this script anyway [yes/no]:")
            if choice == "yes":
                completed = True
                evaluate_stability(
                    xml_file,
                    total_param,
                    param_values,
                    output,
                    temp_output,
                    restart=restart,
                    restart_int=start,
                    end_int=end,
                    saved_data_paths=saved_data_path,
                )
            elif choice == "no":
                completed = True
                print("Script Interrupted.")
            else:
                print("Please enter yes or no.")
        else:
            completed = True
            evaluate_stability(
                xml_file,
                total_param,
                param_values,
                output,
                temp_output,
                restart=restart,
                restart_int=start,
                end_int=end,
                saved_data_paths=saved_data_path,
            )
