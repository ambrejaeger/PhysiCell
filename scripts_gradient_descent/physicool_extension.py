from xml.etree import ElementTree
from scipy import io as sio
import os
import scripts_sensibility_analysis.membrane_growth_functions as ssa

# To avoid having to write the full string for each XML element, I use this simple dict
# The goal is to only specify the name of the parameter (i.g., 'cell_cell_adhesion_strength')
# The get_xml_stem() function will take care of the rest using this dict
CELL_DEFINITIONS_DICT = {
    # TODO: What is the role of the "set relative/absolute" distances?
    # Should these parameters be included?
    'cycle/phase_transition_rates': [
        'rate'
    ],
    'death': [
        "model[@name='apoptosis']/death_rate"
        ],
    'mechanics': [
        'cell_cell_adhesion_strength',
        'cell_cell_repulsion_strength',
        'relative_maximum_adhesion_distance'
    ],
    # TODO: should the options be represented here?
    # For a parameter sweep, I don't think they are relevant because they are Booleans
    # or categorical values, and they change the cell behaviour significantly (2D vs 3D,...)
    'motility': [
        'speed',
        'persistence_time',
        'migration_bias'
    ]
}


def get_cell_xml_stem(key, definition_name='epi_inter'):
    """Returns the XML element name that corresponds to the passed parameter key."""
    for group, params in CELL_DEFINITIONS_DICT.items():
        if key in params:
            return f'cell_definitions/cell_definition[@name="{definition_name}"]/phenotype/{group}/{key}'


def get_me_xml_stem(key, substance='substrate'):
    """Returns the XML element name that corresponds to the passed parameter key."""
    return f'microenvironment_setup/variable[@name={substance}/{key}]'


def update_config_file(params_dict, config_path, cell_type):
    """Updates configuration file with the specified input values."""
    tree = ElementTree.parse(config_path)

    for param, value in params_dict.items():
        param_type, key = param.split('/', maxsplit=1)

        # Expects the structure "cell/key}'
        if param_type == 'cell':
            param_name = get_cell_xml_stem(key,cell_type)

        # Expects the structure "me/{substance}:{key}'
        elif param_type == 'me':
            substance, substance_key = key.split(':')
            param_name = get_me_xml_stem(substance_key, substance)

        #Expects the structure "cell_rule/{rule number}/{parameter}"
        elif param_type == 'cell_rule':
            #find the path to the cell_rule file in the config_file
            cell_rule_file = os.path.join(str(tree.find("cell_rules/rulesets/ruleset/folder").text), str(tree.find("cell_rules/rulesets/ruleset/filename").text))
            #modify the cell rule file
            ssa.modify_csv(cell_rule_file,param_type + '/' + key, value)
            continue
        else:
            param_name = f'custom_variables/{key}'

        tree.find(param_name).text = str(value)

    tree.write(config_path)



def get_cell_data(timestep, folder_name, variables='all'):
    """Returns a dictionary with the cell output data for the selected variables.

    Parameters
    ----------
    timestep : int
        The time point at which the output was recorded
    folder_name: Path
        The path to the folder where the output (.mat) files are stored
    variables : list
        The variables to be extracted from the output files. If variables
        are not defined, all the available outputs will be saved.
    """

    # All possible output variables written by PhysiCell
    data_labels = [
        'ID',
        'position_x', 'position_y', 'position_z',
        'total_volume',
        'cell_type',
        'cycle_model', 'current_phase', 'elapsed_time_in_phase',
        'nuclear_volume', 'cytoplasmic_volume',
        'fluid_fraction', 'calcified_fraction',
        'orientation_x', 'orientation_y', 'orientation_z',
        'polarity',
        'migration_speed',
        'motility_vector_x', 'motility_vector_y', 'motility_vector_z',
        'migration_bias',
        'motility_bias_direction_x', 'motility_bias_direction_y', 'motility_bias_direction_z',
        'persistence_time',
        'motility_reserved'
    ]

    # Create path name
    time_str = str(timestep).zfill(8)
    file_name = 'output{}_cells.mat'.format(time_str)
    path_name = folder_name / file_name

    # Read output file
    cell_data = sio.loadmat(path_name)['cells']

    # Select and save the variables of interest
    variables_indexes = [data_labels.index(var) for var in variables]
    cells = {var: cell_data[index, :]
             for var, index in zip(variables, variables_indexes)}

    return cells