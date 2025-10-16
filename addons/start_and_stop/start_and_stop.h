
#include "../../core/PhysiCell.h"
#include "../../modules/PhysiCell_standard_modules.h" 


using namespace BioFVM; 
using namespace PhysiCell;

void reset_microenv(std::string saved_files_folder);

// Function to save cell data
void save_cell_microenv_data(Cell_Container* cell_container, std::string saved_files_folder);

void reset_cell(double last_cell_cycle_time, std::string saved_files_folder, std::string xml_path);

void reset_global_parameters(Cell_Container* cell_container, std::string saved_files_folder);

void reset_randomness(std::string saved_files_folder);