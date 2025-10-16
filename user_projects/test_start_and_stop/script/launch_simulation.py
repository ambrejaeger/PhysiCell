import sys
import os
sys.path.append('../')
from interface import *

def launch_autostop_simu(two_D : bool, xml_filepath : str = "config/PhysiCell_settings.xml", json_file_path : str = 'script/simulation_parameters/simulation_parameters.json', exe_file : str = "./heterogeneity") -> str:
    root_dir = os.getcwd()
    json_file_path = os.path.join(root_dir, json_file_path)

    # define interface object
    my_interface = Interface(root_dir, json_file_path, two_D)

    # update_parameters
    iteration = 0

    my_interface.update_parameters(iteration, xml_filepath)

    # Execute the simulation with new parameters

    output_folder = my_interface.execute_simulation(iteration, exe_file)
    
    return output_folder


def launch_batch_autostop_simu(num_iter : int, two_D : bool, output_root : str = "./output/", xml_filepath : str = "config/PhysiCell_settings.xml", json_file_path : str = 'script/simulation_parameters/simulation_parameters.json', exe_file : str = "./heterogeneity") -> None:
    root_dir = os.getcwd()
    json_file_path = os.path.join(root_dir, json_file_path)

    # define interface object
    interface = Interface(root_dir, json_file_path, two_D)

    #Execute the first simulation
    iteration = 0
    interface.output_folder = output_root + "output_" + str(iteration)
    interface.parameter_dict["saving_folder"]["value"] = "./start_and_stop_saving_files_" + str(iteration) + "/"

    if not os.path.isdir(interface.output_folder):
            os.mkdir(interface.output_folder)

    interface.execute_simulation(iteration, exe_file)

    for i in range(1, num_iter):
        #Update the name of the output folder
        interface.output_folder = output_root + "output_" + str(i)

        #Update the seed
        interface.parameter_dict["random_seed"]["value"] = i

        #Modify the folder where start_and_stop_saving files are stored 
        interface.parameter_dict["saving_folder"]["value"] = "./start_and_stop_saving_files_" + str(i) + "/"
        interface.parameter_dict["output_folder"]["value"] = output_root + "output_" + str(i)

        if not os.path.isdir(interface.parameter_dict["output_folder"]["value"]):
            os.mkdir(interface.parameter_dict["output_folder"]["value"])

        #update read 
        #Careful here we put iteration which is equal to 0 because we want to generate initial simulation not start from saved file
        interface.update_parameters(iteration, xml_filepath)

        #Run simulation
        #However here we put i because we don't which to remake the files we want to use the same .exe
        interface.execute_simulation(i, "./heterogeneity")

    return None

if __name__ == "__main__":

    os.chdir('.')
    two_D = True
    launch_batch_autostop_simu(2, True)
    #single_simu(two_D)