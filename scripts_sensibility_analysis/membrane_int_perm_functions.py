from membrane_permeability_functions import *
from membrane_integrity_functions import *

def compute_membrane_crossing_time(mat_files, temp_output_folder, tolerance, X_min=-300, X_max=300, save_interval_time=200.0, membrane_type = 2, attracted_type = 4):
    initial_xml_file = os.path.join(temp_output_folder,"initial.xml")
    for i,mat_file in enumerate(mat_files): 
        cells_info = extract_position_type_data(mat_file, initial_xml_file)
        cells_type = cells_info[1] 
        positions = cells_info[0] 

        heights = cell_type_height(cells_info, X_min, X_max, 15, membrane_type)
        attracted = cells_type == attracted_type
        result = cell_below(positions[attracted][0], heights)
        if result:
            print(mat_file)
            return i*save_interval_time
        
    return 0


def evaluate_membrane_int_perm(xml_file, param_treepaths, param_values, tolerance, save_output_folder, temp_output, restart=False, restart_int = 0, end_int = 0): 
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
            np.savetxt(os.path.join(output_folder,"param_values_membrane_int_perm.txt"), param_values)  

        if os.path.exists(xml_file):
            shutil.copy(xml_file, temp_xml_file)
        
        start_file = restart_int
    else:
        if os.path.exists(xml_file) & (not os.path.exists(temp_xml_file)):
            shutil.copy(xml_file, temp_xml_file)

        param_values = np.loadtxt(os.path.join(output_folder,"param_values_membrane_int_perm.txt"))
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

    output_storage_perm = np.zeros(param_values.shape[0])
    output_storage_int = np.zeros(param_values.shape[0])

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
        output_storage_perm[i] = compute_membrane_crossing_time(mat_files, temp_output_folder, tolerance=tolerance, membrane_type = 2, attracted_type = 4)
        output_storage_int[i] = math.fsum(compute_membrane_integrity(mat_files, temp_output_folder, 2, tolerance))
        print("Run ", i, " is done")

        process1 = subprocess.run(
        ["make", "gif", f"OUTPUT={temp_output_folder}"],
        capture_output=True,
        text=True
        )
        
        shutil.copyfile(f"{temp_output_folder}/out.gif", f"{output_folder}/out_{i}.gif")
        with open(output_storage_file, "a") as f:
            f.write(f"{i} {output_storage_perm[i]} {output_storage_int[i]}\n")
        with open(save_output, "a") as f:
            f.write(process0.stdout)
        with open(save_output, "a") as f:
            f.write(process1.stdout)

    #delete temp output at the end of the run
    if os.path.isdir(temp_output_folder):
            shutil.rmtree(temp_output_folder)

    return output_storage_perm, output_storage_int
