from membrane_growth_functions import *
import argparse

param = [["cell_definitions/cell_definition/phenotype/secretion/substrate/secretion_rate", "epi_inter"],
         ["microenvironment_setup/variable/physical_parameter_set/diffusion_coefficient" , "", "", "div_inhib"],
         ["microenvironment_setup/variable/physical_parameter_set/decay_rate", "", "", "div_inhib"]]

param_bounds = [[1.0, 1000.0],
          [500.0, 1500.0],
          [0.05, 0.5]]

cell_rules = []
cell_rules_bounds = []

names = ["_".join(p) for p in param] + ["_".join(r) for r in cell_rules]
num_vars = len(param) + len(cell_rules)
bounds = param_bounds + cell_rules_bounds
xml_file = "./config/PhysiCell_settings.xml"
param_values = define_set_param(num_vars, names, bounds, sample_size=128)

indexes = list(range(1,param_values.shape[0],50))

print(param_values)


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
                evaluate_epi_growth(xml_file, param, param_values, output, temp_output, restart=restart,restart_int=start, end_int=end)
            elif choice == "no":
                completed = True
                print("Script Interrupted.")
            else:
                print("Please enter yes or no.")
        else:
            completed = True
            evaluate_epi_growth(xml_file, param, param_values, output, temp_output, restart=restart,restart_int=start, end_int=end)
