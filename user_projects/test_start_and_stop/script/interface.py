import os
import subprocess
import xml.etree.ElementTree as ET
import json
import sys

from pctk import multicellds
sys.path.append('../')


class Interface:
    def __init__(self, root_dir : str, json_file_path : str, two_D : bool = True):

        # define the directory
        self.root_dir = root_dir
        self.PhysiCell_dir = os.path.join(root_dir)
        self.two_D = two_D

        #define the parameter dict
        self.json_file_path = json_file_path
        with open(json_file_path, 'r') as file:
            self.parameter_dict = json.load(file)
        
        if self.two_D:
            self.parameter_dict['z_min']['value'] = -10
            self.parameter_dict['z_max']['value'] = 10
            self.parameter_dict['use_2D']['value'] = 'true'
        else:
            self.parameter_dict['z_min']['value'] = -500
            self.parameter_dict['z_max']['value'] = 500
            self.parameter_dict['use_2D']['value'] = 'false'

        self.output_folder = os.path.join(self.PhysiCell_dir, 'output')
        self.Start_Stop_folder = os.path.join(self.PhysiCell_dir, 'start_and_stop_saving_files')

    def update_parameters(self, iteration : int, xml_filepath : str) -> None:

        #read_init should always be set to true if it is not the first iteration meaning the simulation is started from saved files
        if iteration == 0:
            self.parameter_dict['read_init']['value'] = 'false'
        else:
            self.parameter_dict['read_init']['value'] = 'true'

        # File path for the physicell settings
        physicell_setting_file = os.path.join(self.root_dir, xml_filepath)
        
        # Upload XML file
        tree = ET.parse(physicell_setting_file)
        root = tree.getroot()
        
        # Extract list of parameters to update
        parameters = list(self.parameter_dict.keys())
        
        for param in parameters:
            tags = self.parameter_dict[param]['path'].split('/')
            new_value = str(self.parameter_dict[param]['value'])

            # CASO SPECIALE: aggiorna tutti i cell_definition/intracellular/...
            if tags[:3] == ['cell_definitions', 'cell_definition', 'phenotype'] and tags[3] == 'intracellular':
                target_tag = tags[-1]  # es: bnd_filename o cfg_filename
                for cd in root.findall('.//cell_definition'):
                    phenotype = cd.find('phenotype')
                    if phenotype is not None:
                        intracellular = phenotype.find('intracellular')
                        if intracellular is not None:
                            target_element = intracellular.find(target_tag)
                            if target_element is not None:
                                target_element.text = new_value
                continue  # salta il resto del loop per questo parametro

            # Altrimenti: gestione standard
            element = root
            for tag in tags:
                if element is not None:
                    # Unsure what this if block does what are the 'variable'
                    if tag == 'variable' and 'name' in self.parameter_dict[param]:
                        found = False
                        for var in element.findall(tag):
                            if var.attrib['name'] == self.parameter_dict[param]['name']:
                                element = var
                                found = True
                                break
                        if not found:
                            element = None
                            break
                    elif tag != 'variable':
                        element = element.find(tag)

            if element is not None:
                element.text = new_value

        # Salva il file alla fine
        tree.write(physicell_setting_file)
        print('Settings updated succesfully!')
        
        return None


    
    def execute_simulation(self, iteration : int, executable_file : str) -> str:
        # Change current working directory
        os.chdir(self.PhysiCell_dir)

        #Check existence of necessary directories 
        if not os.path.isdir(self.Start_Stop_folder):
            try: 
                os.mkdir(self.Start_Stop_folder)
                print(f"Directory '{self.Start_Stop_folder}' created successfully.")

            except FileExistsError:
                print(f"Directory '{self.Start_Stop_folder}' already exists.")
            except PermissionError:
                print(f"Permission denied: Unable to create '{self.Start_Stop_folder}'.")
            except Exception as e:
                print(f"An error occurred: {e}")

        if not os.path.isdir(self.output_folder):
            try: 
                os.mkdir(self.output_folder)
                print(f"Directory '{self.output_folder}' created successfully.")

            except FileExistsError:
                print(f"Directory '{self.output_folder}' already exists.")
            except PermissionError:
                print(f"Permission denied: Unable to create '{self.output_folder}'.")
            except Exception as e:
                print(f"An error occurred: {e}")

        #If it is the first execution of the simulation then make the project           
        if iteration == 0:

            # Recreate output folder
            first_make_command = ["make", 'data-cleanup']
            subprocess.run(first_make_command, check=True)

            reset_make_command = ["make", 'reset']
            subprocess.run(reset_make_command, check=True)
            
            clean_make_command = ["make", 'clean']
            subprocess.run(clean_make_command, check=True)

            make_command = ["make", "load", "PROJ=test_start_and_stop"]
            subprocess.run(make_command, check=True)
        
            make_command = ["make"]
            subprocess.run(make_command, check=True)
        
        # Run the simulation
        execute_command = ["./" + executable_file]
        subprocess.run(execute_command, check=True) 
        
        return self.output_folder

    



