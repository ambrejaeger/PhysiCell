from sensibility_analysis_functions import *
import argparse

param = [["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_basal", "epi_basal"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_basal", "epi_inter"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "epi_basal", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "epi_basal"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "membrane"],
    ["cell_definitions/cell_definition/phenotype/mechanics/cell_adhesion_affinities/cell_adhesion_affinity", "membrane", "conjonctif"]
    ]

bounds = [[0.0, 1.0] * len(param)]
names = ["_".join(p) for p in param]
num_vars = len(param)
xml_file = "./config/PhysiCell_settings.xml"
param_values = define_set_param(num_vars, names, bounds)
tolerance = 5.0

indexes = list(range(1,param_values.shape[0],50))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that runs evaluate_membrane_integrity2"
    )
    parser.add_argument("--output", required=True, type=str)
    parser.add_argument("--restart", required=False, type=str, choices=["0", "1", "true", "false"])
    parser.add_argument("--start", required=False, type=int)
    parser.add_argument("--end", required=False, type=int)

    args = parser.parse_args()

    output = args.output
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
        
    evaluate_membrane_integrity2(xml_file, param, param_values, tolerance, output,restart=restart,restart_int=start, end_int=end)