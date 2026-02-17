from sensibility_analysis_functions import *
from membrane_growth_functions import *

def compute_median_apoptosis_position(label_file, mat_files):
    dead_y_pos = {} 
    
    for file in mat_files:
        df = get_cells_data(label_file, file)

        
        dead_mask = df.loc['dead'] == 1
        dead_cells_df = df.loc[:, dead_mask]

        for cell_id in dead_cells_df.columns:
            if cell_id not in dead_y_pos:
                dead_y_pos[cell_id] = dead_cells_df.at["position_2", cell_id]

    df = get_cells_data(label_file, mat_files[-1])           
    mask_1 = df.loc['cell_type'] == 1
    median_type_1 = df.loc['position_2', mask_1].median()

    mask_0 = df.loc['cell_type'] == 0
    median_type_0 = df.loc['position_2', mask_0].median()

    return len(dead_y_pos), np.median(list(dead_y_pos.values())), median_type_1, median_type_0

if __name__ == "__main__":

    param = [
    [
        "cell_definitions/cell_definition/phenotype/death/model/death_rate",
        "epi_inter",
        "",
        "",
        "",
    ]
    ]
    #Parameters that were initialized differently when generating the pre-epithelium
    saved_data_path = [["Secretion/Secretion 1/Secretion_Rate", "type", 1], [], [], ["Death/Model 0/Rate", "type", 1]]

    param_bounds = [[0.000001, 0.001]]

    #We wish to look at the effect of the hill power and the half max value of th third rule [half max value, hill power]
    cell_rules = [['cell_rule',"2", "5"], ['cell_rule',"2", "6"]]
    cell_rules_bounds = [[0.0, 6.0], [1.0, 5.0]]


    names = ["_".join(p) for p in param] + ["_".join(str(r)) for r in cell_rules]
    num_vars = len(param) + len(cell_rules)
    bounds = param_bounds + cell_rules_bounds
    xml_file = "./config/PhysiCell_settings.xml"
    param_values = define_set_param(num_vars, names, bounds, sample_size=64)

    total_param = param + cell_rules
    
    result_file = "./output_sensibility_analysis/membrane_stability/epi_growth.txt"
    process_file(result_file, result_file)
    param_names_file = "./output_stable/param_names.txt"
    param_values_file = "./output_stable/param_values.txt"

    #Combine file with header
    output_path = "./output_sensibility_analysis/membrane_stability/result_summary_growth.txt"
    #combine_files_with_header_2(result_file, param_values_file, output_path, names, 'growth_rate', 'epi_size')


    Si = analyze_sobol(result_file, param_names_file, bounds, column=1)
    save_dataframes_to_txt([Si.to_df()[0], Si.to_df()[1], Si.to_df()[2]], "./output_sensibility_analysis/membrane_stability/data_output_growth.txt")
    #Si.plot()
    #plt.show()

    #Plotting scatter plot groupped variable 
    X, Y = extract_X_Y(output_path, names[3], Y_column='growth_rate')
    plot_scatter_sets(Y, X, set_names="", xlabel="Death rate", title ="Epithelium growth rate in function of death rate", ylabel="Growth rate")
"""x
    param = [
        [
            "cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity",
            "epi_inter",
            "epi_inter",
            "",
            "",
        ]
    ]

    saved_data_path = [["Mechanics/cell_adhesion_affinities[1]", "type", 1]]

    param_bounds = [[0.0, 1]]

    cell_rules = [['cell_rule',"2", "5"]]
    cell_rules_bounds = [[0.01, 4]]

    names = ["_".join(p) for p in param] + ["_".join(str(r)) for r in cell_rules]
    num_vars = len(param) + len(cell_rules)
    bounds = param_bounds + cell_rules_bounds
    xml_file = "./config/PhysiCell_settings.xml"
    param_values = define_set_param(num_vars, names, bounds, sample_size=64)

    total_param = param + cell_rules

    result_file = "./output_sensibility_analysis/membrane_apoptosis/epi_growth.txt"
    process_file(result_file, result_file)
    param_names_file = "./output_apop/param_names.txt"
    param_values_file = "./output_apop/param_values.txt"

    #Combine file with header
    output_path = "./output_sensibility_analysis/membrane_apoptosis/result_summary_growth.txt"
    #combine_files_with_header_2(result_file, param_values_file, output_path, names, 'nbr_apop', 'apop_pos', 'epi_size')


    #Si = analyze_sobol(result_file, param_names_file, bounds, column=1)
    #save_dataframes_to_txt([Si.to_df()[0], Si.to_df()[1], Si.to_df()[2]], "./output_sensibility_analysis/membrane_apoptosis/data_output_growth.txt")
    #Si.plot()
    #plt.show()

    #Plotting scatter plot groupped variable 
    #X, Y = extract_X_Y(output_path, names[0], Y_column='epi_size')
    #plot_scatter_sets(Y, X, set_names="", xlabel="adhesion", title ="Epithelium growth rate in function of adhesion", ylabel="Epithelium size")
    """