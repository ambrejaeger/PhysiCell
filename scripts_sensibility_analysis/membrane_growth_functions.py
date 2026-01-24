from sensibility_analysis_functions import *

#Write a function to compute the thickness of the epi inter, and its growth rate
def compute_epi_thickness(mat_file, xml_file):
    positions,types,cells_ID =  extract_position_type_data(mat_file, xml_file)
    
    mask_epi_inter_cells = (types == 1.0)
    positions_epi_inter = positions[mask_epi_inter_cells]
    positions_y_epi_inter = np.sort(positions_epi_inter[:, 1])
    
    min = np.average(positions_y_epi_inter[:10])
    max = np.average(positions_y_epi_inter[-10:])

    return abs(max - min)

#Write a function to test if the epithelium is stable in size
def compute_epi_stability(mat_files, xml_file, sort = False):
    if (sort):
        #Sorting only works if your files are named acccording to PhysiCell standard outputs
        mat_files.sort(key=lambda x: int(x.split('_')[0].split('output')[-1]))
        print(mat_files)

    sizes = []
    for file in mat_files:
        sizes.append(compute_epi_thickness(file, xml_file))

    growth_rates = []
    for i in range(1, len(sizes)):
        growth_rates.append(sizes[i] - sizes[i-1])
    return np.average(growth_rates)

#Write a function to track the number of the cells type
def compute_number_cells(mat_file, xml_file, tracked_types):
    positions,types,cells_ID =  extract_position_type_data(mat_file, xml_file)
    types_dict = {}

    for t in tracked_types:
        mask = (types == t)
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
#Write a function to evaluate the membrane thickness


def modify_csv(file_path, cell_rule, value):
    # Validate cell_rule format
    if len(cell_rule) != 3 or cell_rule[0] != 'cell_rules':
        raise ValueError("cell_rule must be in format: ['cell_rule', row_index, col_index]")
    
    # Get row and column indices
    row_idx = cell_rule[1]
    col_idx = cell_rule[2]
    
    # Read all lines from the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Check if row index is valid
    if row_idx < 0 or row_idx >= len(lines):
        raise IndexError(f"Row index {cell_rule[1]} is out of range. File has {len(lines)} rows.")
    
    # Split the line by commas and modify the specific column
    line_parts = lines[row_idx].strip().split(',')
    
    # Check if column index is valid
    if col_idx < 0 or col_idx >= len(line_parts):
        raise IndexError(f"Column index {cell_rule[2]} is out of range. Row {cell_rule[1]} has {len(line_parts)} columns.")
    
    # Modify the specific cell
    line_parts[col_idx] = str(value)
    
    # Reconstruct the line
    lines[row_idx] = ','.join(line_parts) + '\n'
    
    # Create the output file path
    folder = os.path.dirname(file_path)
    output_path = os.path.join(folder, 'temp_cell_rules.csv')
    
    # Write all lines to the new file
    with open(output_path, 'w') as file:
        file.writelines(lines)
    
    print(f"Modified copy saved to: {output_path}")
    print(f"Modified cell at row {cell_rule[1]}, column {cell_rule[2]} to value: {value}")
    
    return output_path


def evaluate_epi_growth(xml_file, param_treepaths, param_values, save_output_folder, temp_output, restart=False, restart_int = 0, end_int = 0): 
    output_folder = os.path.join(os.getcwd(), save_output_folder)
    save_output = os.path.join(output_folder, "run_output.txt")
    output_storage_file = os.path.join(output_folder,"epi_growth.txt")
    output_storage_cell_pop = os.path.join(output_folder,"epi_cell_pop.txt")
    temp_output_folder = os.path.join(os.getcwd(), temp_output)
    temp_xml_file = os.path.join(temp_output_folder, os.path.basename(xml_file))
    start_file = 0

    #Identify the cell_rules file if defined
    tree = ET.parse(xml_file)
    root = tree.getroot()
    filepath = ''
    if root.find('./cell_rules/rulesets/ruleset/folder') != None:
        folder = root.find('./cell_rules/rulesets/ruleset/folder').text 
        file = root.find('./cell_rules/rulesets/ruleset/filename').text
        print(file)
        print(folder)
        filepath = os.path.join(folder,file)
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
            with open(os.path.join(output_folder, "param_names.txt"), 'w') as fp:
                for item in param_treepaths:
                    fp.write("%s\n" % item)
            np.savetxt(os.path.join(output_folder,"param_values.txt"), param_values)  

        if os.path.exists(xml_file):
            shutil.copy(xml_file, temp_xml_file)
        
        start_file = restart_int
    else:
        if os.path.exists(xml_file) & (not os.path.exists(temp_xml_file)):
            shutil.copy(xml_file, temp_xml_file)

        param_values = np.loadtxt(os.path.join(output_folder,"param_values.txt"))
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
    if not os.path.isfile(output_storage_cell_pop):
        with open(output_storage_cell_pop, 'w') as f:
            pass
    

    output_growth_rates = []
    output_epi_sizes = []
    output_cell_pop = []
    #Changing output folder
    modify_xml(temp_xml_file, "save/folder", temp_output)

    for i in range(end_file - start_file): 
        for j, val in enumerate(param_values[i + start_file,:]):
            name_cell_def = ''
            name_interact_cell_def = ''
            variable_name = ''
            substrate = ''

            if len(param_treepaths[j]) > 1:
                name_cell_def = param_treepaths[j][1]
            if len(param_treepaths[j]) > 2:
                name_interact_cell_def = param_treepaths[j][2]
            if len(param_treepaths[j]) > 3:
                variable_name = param_treepaths[j][3]
            if len(param_treepaths[j]) > 4:
                substrate = param_treepaths[j][4]


            if param_treepaths[j][0] == "cell_rules":
                modify_csv(filepath,param_treepaths[j], val)
                modify_xml(temp_xml_file, './cell_rules/rulesets/ruleset/filename', 'temp_cell_rules.csv')
            else:
                modify_xml(temp_xml_file, param_treepaths[j][0], 
                    val, name_cell_def=name_cell_def, 
                    name_interact_cell_def=name_interact_cell_def,
                    variable_name = variable_name,
                    substrate = substrate) 

        #Running simulation 
        process0 = subprocess.run(
                        ["./test_death", temp_xml_file],
                        capture_output=True, 
                        text=True)

        #Outputs analysis
        label_file = os.path.join(temp_output_folder, 'initial.xml')
        mat_files = get_output_files(temp_output_folder)
        #Average growth rate of epi_inter layer over the run
        output_growth_rates.append(np.average(compute_epi_stability(mat_files, label_file)))
        #Size of the epi_inter layer at the last time step
        output_epi_sizes.append(compute_epi_thickness(mat_files[-1], label_file))
        #Cells population size
        output_cell_pop.append(compute_number_cells_over_time(mat_files, label_file, [0, 1]))
        print("Run ", i, " completed")

        process1 = subprocess.run(
                        ["make", "gif", f"OUTPUT={temp_output_folder}"],
                        capture_output=True,
                        text=True)
        
        shutil.copyfile(f"{temp_output_folder}/out.gif", f"{output_folder}/out_{i + start_file}.gif")
        with open(output_storage_file, "a") as f:
            f.write(f"{i + start_file} {output_growth_rates[i]} {output_epi_sizes[i]}\n")
        with open(output_storage_cell_pop, "a") as f:
            f.write(f"{i + start_file} {output_cell_pop[i]}\n")
        with open(save_output, "a") as f:
            f.write(process0.stdout)
        with open(save_output, "a") as f:
            f.write(process1.stdout)

    #delete temp output at the end of the run
    if os.path.isdir(temp_output_folder):
        shutil.rmtree(temp_output_folder)

    return output_growth_rates, output_epi_sizes

def main():
    #mat_file = "./output/output00000110_cells.mat"
    #mat_files = ["./output/output00000109_cells.mat", "./output/output00000108_cells.mat", "./output/output00000110_cells.mat", "./output/output00000111_cells.mat", "./output/output00000112_cells.mat"]
    xml_file = "./config/PhysiCell_settings.xml"
    #print(compute_epi_thickness(mat_file, xml_file))
    #print(compute_epi_stability(mat_files, xml_file, sort = True))
    #compute_number_cells_over_time(mat_files, xml_file, [0, 1])
    param = [['cell_rules',1,5]]
    value = 4.0
    #modify_xml(xml_file, path, value, name_cell_def=param[0][1], name_interact_cell_def="", variable_name="", substrate=param[0][4])
    filepath = './config/cell_rules.csv'

    modify_csv(filepath, param[0],0.5)

if __name__ == "__main__":
    main()