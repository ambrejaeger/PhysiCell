from physicool import optimization as opt
from physicool_v2.updaters import CellUpdater, update_all
from physicool_v2.processing import get_final_y_position

# Compiles the project and creates a black box object for it
# opt.compile_project()    
black_box = opt.PhysiCellBlackBox(project_name="test_death")

# Define the updater we want to use (change motility data)
new_values = {"speed": 2.0, "migration_bias": 0.9}
updater = CellUpdater(updater_function=update_all,
                     config_path="config/PhysiCell_settings.xml", cell_definition_name="epi_inter")

# Assign the updater and processor to the black box
black_box.updater = updater
black_box.processor = get_final_y_position
black_box.version = "1.9.1"

# Run the model with the target values (speed=2.0 bias=0.9)
opt.clean_tmp_files()
target_data = black_box.run(new_values)
