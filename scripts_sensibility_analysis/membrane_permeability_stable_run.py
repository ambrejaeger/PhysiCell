from sensibility_analysis_functions import *
import argparse

from membrane_int_perm_functions import *

#Parameters to test are speed, migration bias, and cell volume for attracted cell and adhesion affinity between membrane cells, ruser_projects/test_permeability/script user_projects/test_start_and_stop/script/membrane_integrity_run.py user_projects/test_start_and_stop/script/sensibility_analysis_functions.pyepulsion and adhesion between membrane and attracted cell

param = [["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "attracted"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "attracted", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "conjonctif"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "conjonctif", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "conjonctif", "conjonctif"],
    ["cell_definitions/cell_definition/phenotype/motility/speed", "attracted"]
    ]
""",
    ["cell_definitions/cell_definition/phenotype/volume/total", "attracted"],
    ["cell_definitions/cell_definition/phenotype/mechanics/attachment_rate", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/attachment_rate", "conjonctif"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_cell_repulsion_strength", "attracted"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_cell_repulsion_strength", "membrane"]
    
    ,
          [3000, 5500],
          [0.0, 1.0],
          [0.0, 1.0],
          [10.0, 100.0],
          [10.0, 100.0]"""

bounds = [[0.0, 1.0],
          [0.0, 1.0],
          [0.0, 1.0],
          [0.0, 1.0],
          [0.0, 1.0],
          [0.0, 1.0],
          [1.0, 5.0]
          ]
names = ["_".join(p) for p in param]
num_vars = len(param)
xml_file = "./config/PhysiCell_settings.xml"
groups = ['Group_m1', 'Group_m2', 'Group_m2', 'Group_m3', 'Group_m3', 'Group_c1', 'Group_sattr', 'Group_vattr', 'Group_m_att', 'Group_c_att', 'Group_at_rep', 'Group_m_rep']
param_values = define_set_param(num_vars, names, bounds, groups = groups, sample_size = 256)
tolerance = 5.0

indexes = list(range(1,param_values.shape[0],50))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that runs evaluate_membrane_integrity2"
    )
    parser.add_argument("--output", required=True, type=str)
    parser.add_argument("--temp_output", required=True, type=str)
    parser.add_argument("--restart", required=False, type=str, choices=["0", "1", "true", "false"])
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
            print(temp_output, " folder already exist in your working directory, by running this script you might overwrite an ongoing simulation.")
            choice = input("Do you want to run this script anyway [yes/no]:")
            if choice == "yes":
                completed = True
                evaluate_membrane_int_perm(xml_file, param, param_values, tolerance, output, temp_output, restart=restart,restart_int=start, end_int=end)
            elif choice == "no":
                completed = True
                print("Script Interrupted.")
            else:
                print("Please enter yes or no.")
        else:
            completed = True
            evaluate_membrane_int_perm(xml_file, param, param_values, tolerance, output, temp_output, restart=restart,restart_int=start, end_int=end)
