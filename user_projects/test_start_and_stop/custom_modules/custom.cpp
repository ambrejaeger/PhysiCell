/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include <algorithm>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <cstdlib>
#include <cmath>
#include <queue>
#include <functional>
#include <sstream>
#include <vector>
#include <string>
#include <deque>
#include <unordered_set>
#include <filesystem>
#include <sys/stat.h>

#include "./custom.h"
#include "../modules/PhysiCell_geometry.h"

void create_cell_types( void )
{
	// Setting the random seed in load_xml_file as of PhysiCell 1.14.0
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	
	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	cell_defaults.functions.division_orientation = custom_division_orientation;

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 
	if (parameters.bools("outer_box"))
	{
		std::string token;
		std::stringstream ss(parameters.strings("outer_box_bounds"));
		std::vector<float> result;
		while (std::getline( ss, token, ',')) {
        	result.push_back(std::stof(token));
		}
		std::cout << result[0] << "," << result[1] << std::endl;
		microenvironment.mesh.outer_bounding_box[0] = result[0];
		microenvironment.mesh.outer_bounding_box[1] = result[1];
		microenvironment.mesh.outer_bounding_box[2] = result[2];
		microenvironment.mesh.outer_bounding_box[3] = result[3];
		microenvironment.mesh.outer_bounding_box[4] = result[4];
		microenvironment.mesh.outer_bounding_box[5] = result[5];
	}
	if (evaluate_start_stop_parameters() != -1) {
		std::cout << "Evaluate start stop parameters completed" << std::endl;
	}
		
	
	return; 
}

void setup_tissue( void )
{
	std::vector<std::vector<double>> positions;
	//Check first if cells are initialized as to end of other simulation
	if ( parameters.bools("read_init") )
	{
		std::string csv_fname = parameters.strings("init_cells_filename");
		positions = read_cells_positions(csv_fname, '\t', true);
		if (positions.empty()) 
		{
			std::cout << "Unproper initialization from " << csv_fname << std::endl;
			return;
		}
		//Creating cells at position from init.tsv or random positions
		Cell* pC; 
		for (int i = 0; i < positions.size(); i++)
		{

			pC = create_cell(get_cell_definition("default"));
			pC->assign_position(positions[i]);
		}

		std::cout << positions.size() << " cells positioned from previous simulation" << std::endl;
		return; 
	}
	else
	{
		//Check if cells have to be initialized by csv
		bool loaded = load_cells_from_pugixml();
		if (loaded) 
		{ 
			std::cout << "Cells loaded from .csv" << std::endl;
			return; 
		}
		else
		{
			double Xmin = microenvironment.mesh.bounding_box[0]; 
			double Ymin = microenvironment.mesh.bounding_box[1]; 
			double Zmin = microenvironment.mesh.bounding_box[2]; 

			double Xmax = microenvironment.mesh.bounding_box[3]; 
			double Ymax = microenvironment.mesh.bounding_box[4]; 
			double Zmax = microenvironment.mesh.bounding_box[5]; 
		
			if( default_microenvironment_options.simulate_2D == true )
			{
				Zmin = 0.0; 
				Zmax = 0.0; 
			}
		
			double Xrange = Xmax - Xmin; 
			double Yrange = Ymax - Ymin; 
			double Zrange = Zmax - Zmin; 
		
			// create some of each type of cell 
			Cell* pC;

			for( int k=0; k < cell_definitions_by_index.size() ; k++ )
			{
				Cell_Definition* pCD = cell_definitions_by_index[k];  
				for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
				{
					std::vector<double> position = {0,0,0}; 
					position[0] = Xmin + UniformRandom()*Xrange; 
					position[1] = Ymin + UniformRandom()*Yrange; 
					position[2] = Zmin + UniformRandom()*Zrange; 
					
					pC = create_cell( *pCD ); 
					pC->assign_position( position );
					
				}
			}
		}
	}
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

std::vector<std::string> heterogeneity_coloring_function( Cell* pCell )
{
	double p = get_single_signal( pCell, "custom:oncoprotein"); 
	
	static double p_min = parameters.doubles( "oncoprotein_min" ); 
	static double p_max = parameters.doubles( "oncoprotein_max" ); 
	
	// immune are black
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ return output; } 
	
	// live cells are green, but shaded by oncoprotein value 
	if( pCell->phenotype.death.dead == false )
	{
		int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (p-p_min) * 255.0 ); 
		char szTempString [128];
		sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein );
		output[0].assign( szTempString );
		output[1].assign( szTempString );

		sprintf( szTempString , "rgb(%u,%u,%u)", (int)round(output[0][0]/p_max) , (int)round(output[0][1]/p_max) , (int)round(output[0][2]/p_max) );
		output[2].assign( szTempString );
		
		return output; 
	}

	// if not, dead colors 
	
	if( get_single_signal( pCell, "apoptotic") > 0.5 )
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	if( get_single_signal(pCell, "necrotic") > 0.5 )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
}

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt )
{
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// if cell is dead, don't bother with future phenotype changes. 
	if( get_single_signal( pCell, "dead") > 0.5 )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

	// multiply proliferation rate by the oncoprotein 

	double cycle_rate = get_single_behavior( pCell, "cycle entry"); 
	cycle_rate *= get_single_signal( pCell , "custom:oncoprotein"); 
	set_single_behavior( pCell, "cycle entry" , cycle_rate ); 
	
	return; 
}


/*******************************************/
/*  FUNCTIONS FOR SPHEROID INITIALIZATION  */
/*******************************************/

std::vector<std::vector<double>> read_cells_positions(std::string filename, char delimiter, bool header)
{
	// File pointer
	std::fstream fin;
	std::vector<std::vector<double>> positions;

	// Open an existing file
	fin.open(filename, std::ios::in);

	// Read the Data from the file
	// as String Vector
	std::vector<std::string> row;
	std::string line, word;

	if (header)
	{ getline(fin, line); }

	do
	{
		row.clear();

		// read an entire row and
		// store it in a string variable 'line'
		getline(fin, line);
		if (line.empty()) {continue;}

		// used for breaking words
		std::stringstream s(line);

		while (getline(s, word, delimiter))
		{ 
			row.push_back(word);
		}
		try
		{
			if (row.size() == 3)
			{ 
				std::vector<double> tempPoint(3,0.0);
				tempPoint[0]= std::stof(row[0]);
				tempPoint[1]= std::stof(row[1]);
				tempPoint[2]= std::stof(row[2]);

				positions.push_back(tempPoint);
			}
			else
			{
				throw std::runtime_error("");
			}
		}
		catch(const std::exception& e)
		{
			std::cerr << "Improper formatting of " << filename << " - possible empty lines or invalid data. Error: " << e.what() << std::endl;
			positions.clear();
			return positions;
		}

	} while (!fin.eof());
	
	return positions;
}


std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);
	
	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);
	
	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;
				
				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
			
		}
	}
	return cells;
	
}


std::vector<std::vector<double>> create_cell_disc_positions(double cell_radius, double disc_radius)
{	 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	double x = 0.0; 
	double y = 0.0; 
	double x_outer = 0.0;

	std::vector<std::vector<double>> positions;
	std::vector<double> tempPoint(3,0.0);
	
	int n = 0; 
	while( y < disc_radius )
	{
		x = 0.0; 
		if( n % 2 == 1 )
		{ x = 0.5 * cell_spacing; }
		x_outer = sqrt( disc_radius*disc_radius - y*y ); 
		
		while( x < x_outer )
		{
			tempPoint[0]= x; tempPoint[1]= y;	tempPoint[2]= 0.0;
			positions.push_back(tempPoint);			
			if( fabs( y ) > 0.01 )
			{
				tempPoint[0]= x; tempPoint[1]= -y;	tempPoint[2]= 0.0;
				positions.push_back(tempPoint);
			}
			if( fabs( x ) > 0.01 )
			{ 
				tempPoint[0]= -x; tempPoint[1]= y;	tempPoint[2]= 0.0;
				positions.push_back(tempPoint);
				if( fabs( y ) > 0.01 )
				{
					tempPoint[0]= -x; tempPoint[1]= -y;	tempPoint[2]= 0.0;
					positions.push_back(tempPoint);
				}
			}
			x += cell_spacing; 
		}		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
	return positions;
}

void inject_density_sphere(int density_index, double concentration, double membrane_lenght)
{
	// Inject given concentration on the extremities only
	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		auto current_voxel = microenvironment.voxels(n);
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};

		if ((membrane_lenght - norm(cent)) <= 0)
			microenvironment.density_vector(n)[density_index] = concentration;
	}
}

void remove_density(int density_index)
{
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
		microenvironment.density_vector(n)[density_index] = 0;
}


double total_live_cell_count()
{
        double out = 0.0;

        for( int i=0; i < (*all_cells).size() ; i++ )
        {
                if( (*all_cells)[i]->phenotype.death.dead == false && (*all_cells)[i]->type == 0 )
                { out += 1.0; }
        }

        return out;
}

double total_dead_cell_count()
{
        double out = 0.0;

        for( int i=0; i < (*all_cells).size() ; i++ )
        {
                if( (*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 0 )
                { out += 1.0; }
        }

        return out;
}


/****************************************/
/* START AND STOP FUNCTIONS DEFINITIONS */
/****************************************/

using namespace std;

vector<double> vector_alives;

std::unordered_map<std::string, bool> auto_stop_param = {{"start_stop", false}, {"read_init", false},{"saving_folder", false}, {"init_cells_filename", false}, {"auto_stop", false}, {"auto_stop_alive", false}, {"auto_stop_epi_stable", false}, {"auto_stop_epi_size", false}};

int evaluate_start_stop_parameters() {
	int result = 1;
	//auto stop parameters should be boolean or string in the user_parameters section of your .xml config files
	std::cout << "Auto stop user parameters evaluation: " << std::endl;
	
	if (parameters.bools.size() > 0) {
		for ( auto &p : auto_stop_param ) {
			if ( parameters.bools.find_index(p.first) != -1 ) {
				p.second = parameters.bools(p.first); 
				std::cout << p.first << " is " << p.second << std::endl;
			}
		}
	}

	//Check that init file and saving folder exists
	if (auto_stop_param["start_stop"]) {
		if (parameters.strings.size() > 0) {
			if ( auto_stop_param["read_init"] ) {
				std::cout << "This is running 1" << std::endl;
				if ( parameters.strings.find_index("init_cells_filename") != -1 ) {
					std::cout << "This is running 2" << std::endl;
					auto_stop_param["init_cells_filename"] = true;
					const char *file = parameters.strings("init_cells_filename").c_str();
					struct stat sb;
					if (stat(file, &sb) == 0 && !(sb.st_mode & S_IFDIR)) {
						std::cout << "The init file: " << file << " exists" << std::endl;
					}
					else {
						std::cout << "The file at the path " << parameters.strings("init_cells_filename") << " does not exist" << std::endl;
						result = -1;
					}
				}
				else {
					std::cout << "Path to the init file not indicated in <user_parameters> with the tag <init_cells_filename> in the .xml config file" << std::endl;
					result = -1;
				}
			}
			if ( parameters.strings.find_index("saving_folder") != -1 ) {
				auto_stop_param["saving_folder"] = true;
				const char *dir = parameters.strings("saving_folder").c_str();
				struct stat sb;
				if (stat(dir, &sb) == 0 ) {
					std::cout << "The directory " << dir << " exists" << std::endl;
				}
				else {
					std::cout << "The directory at the path " << parameters.strings("saving_folder") << " does not exist" << std::endl;
					result = -1;
				}
			}
			else {
				std::cout << "Path to the saving folder is not indicated in <user_parameters> with the tag <saving_folder> in the .xml config file" << std::endl;
				result = -1;
			}
		}
	}
	return result; 
}

bool auto_stop_epi_size(vector<double> vector_epi_pos, double epi_max_size) {
    //concatenate the number of alive cells to the vector
	double epi_average_size = accumulate(vector_epi_pos.begin(), vector_epi_pos.end(), 0);
	epi_average_size /= vector_epi_pos.size();

	std::cout << "Epithelium thickness: " << epi_average_size << std::endl;

	bool condition = false;
	// compute the derivative only for the last three steps
	if (epi_average_size >= epi_max_size) {
		condition = true;
	} else {
		condition = false;
	}

	bool stop;

    if (condition) {
        stop = true;
    } else {
        stop = false;
    }
    return stop;
}


bool auto_stop_epi_stable(vector<double> vector_epi_pos, deque<double> &deque_epi_average_size, double steps, double tolerance)
{
	double epi_average_size = accumulate(vector_epi_pos.begin(), vector_epi_pos.end(), 0);
	epi_average_size /= vector_epi_pos.size();
	deque_epi_average_size.push_back(epi_average_size);

	bool stable = true;
	if ( deque_epi_average_size.size() == steps )
	{
		for (int i = 1; i < deque_epi_average_size.size(); i++)
		{
			if (abs(deque_epi_average_size[i] - deque_epi_average_size[i-1]) > tolerance )
			{
				stable = false;
			}
		}
	}
	else
	{
		stable = false;
	}

	bool stop;

    if (stable) {
        stop = true;
		std::cout << "auto stop stable epithelium size condition activated, simulation interrupted." << std::endl;
    } else {
        stop = false;
		if (deque_epi_average_size.size() >= steps)
		{
			deque_epi_average_size.pop_front();
		}
    }
    return stop;

}

bool auto_stop_alive(int alive_cells) {
    //concatenate the number of alive cells to the vector
	vector_alives.push_back(alive_cells);
	std::cout << "Steps: " << vector_alives.size() << std::endl;

	bool condition = false;
	// check the number of elements inside the vector to decide if process it and compute the derivative
	if (vector_alives.size() >= 8) {
		std::vector<double> derivative;

		// compute the derivative only for the last three steps
		for (size_t i = vector_alives.size() - 4; i < vector_alives.size(); ++i) {
			double slope = vector_alives[i] - vector_alives[i - 1];
			derivative.push_back(slope);
		}

		condition = true;
		for (double slope : derivative) {
			// if the slope is less than or equal to zero, set condition to false
			if (slope > 100) {
				condition = false;
				break;
			}
		}

	} else {
		condition = false;
	}

	bool stop;

    if (condition) {
        stop = true;
    } else {
        stop = false;
    }
    return stop;
}


bool auto_stop() {

	bool condition = false;
	bool stop;
	// implement here your condition to stop the simulation

    if (condition) {
        stop = true;
    } else {
        stop = false;
    }
    return stop;
}

/*****************************************/
/*  FUNCTIONS TO CREATE PRE-EPITHELIUM  */
/*****************************************/
void create_pre_epithelium( int argc, char* argv[] ) {
		
		std::string xml_path_str = argv[1];	
		std::cout << "This runs 2" << std::endl;
		load_PhysiCell_config_file( xml_path_str);
		std::cout << "This runs 3" << std::endl;
		std::string time_units = "min"; 

		/* Microenvironment setup */ 
	
		setup_microenvironment();
		double mechanics_voxel_size = 30; 
		Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

		/* Users typically start modifying here. START USERMODS */ 
		create_cell_types();

		std::cout << "This runs 4" << std::endl;

		std::vector<std::vector<double>> positions;
	
		//Create nbr of .csv for epithelium initialization
		int nbr = 1;
		std::string func = "default";

		if ( argc > 2 ) { nbr = std::stoi(argv[2]); }
		if ( argc > 3 ) { func = argv[3]; }
		create_epithelium_csv(nbr, func);

		return ;
}

void create_epithelium_csv(int nbr, std::string func ) {
	//add posibility to differentiate between 2D and 3D
	for (int i = 0; i < nbr; i++) {
		if (func == "custom"){
		
		}
		else { 
			std::string filename = "cells_" + std::to_string(i) + ".csv";
			position_epithelium_cells();
			save_cells_csv(filename);
		}
	}
	return;
}

void position_epithelium_cells() {
	//Function for 2D pre-epithelium creation

	std::vector<std::string> epi_cell_types = {"epi_basal", "conjonctif", "membrane"};
	
	try{
		for (auto type : epi_cell_types)
		{
			if (cell_definition_indices_by_name.find(type) == cell_definition_indices_by_name.end())
			{
				throw std::runtime_error("Missing cell type: " + type);
			}
		}
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}

	Cell_Definition *pCD_c = find_cell_definition("conjonctif");
	Cell_Definition *pCD_e = find_cell_definition("epi_basal");
	Cell_Definition *pCD_m = find_cell_definition("membrane");

	double conjonctive_layer_thickness = 100;
	double epi_basal_layer_thickness = pCD_e->phenotype.geometry.radius * 2;
	double membrane_layer_thickness = pCD_m->phenotype.geometry.radius * 4;

	//Get size of the microenvironment
	double xmin = default_microenvironment_options.X_range[0];
	double xmax = default_microenvironment_options.X_range[1];
	double ymin = default_microenvironment_options.Y_range[0];
	double ymax = default_microenvironment_options.Y_range[1];


	//Creating conjonctive layer
	std::vector<double> bounds_c = {xmin,ymin,xmax,ymin + conjonctive_layer_thickness};
	random_fill_rectangle(bounds_c, pCD_c);

	//Creating membrane layer
	std::vector<double> bounds_m = {xmin, ymin + conjonctive_layer_thickness, 0, xmax, ymin + conjonctive_layer_thickness + membrane_layer_thickness, 0};
	fill_rectangle(bounds_m, pCD_m);

	//Creating epi_basal layer
	std::vector<double> bounds_e = {xmin,ymin + conjonctive_layer_thickness + membrane_layer_thickness - (pCD_e->phenotype.geometry.radius * 0.5), 0, xmax, ymin + conjonctive_layer_thickness + membrane_layer_thickness + epi_basal_layer_thickness, 0};
	fill_rectangle(bounds_e, pCD_e);
	
	
	return;
}

void random_fill_rectangle (BioFVM::gradient bounds, PhysiCell::Cell_Definition *pCD, double confluence ) {
	//confluence is the proportion of the microenvironment / part of the microenvironment surface occupied by cells

	double cell_radius = pCD->phenotype.geometry.radius;
	double cell_surface = cell_radius * cell_radius * M_PI;
		
	double Xrange = bounds[2] - bounds[0]; 
	double Yrange = bounds[3] - bounds[1]; 

	int number_of_cells;
	double rectangle_surface = Xrange * Yrange;

	number_of_cells = std::round(rectangle_surface / cell_surface * confluence);

	// create some of each type of cell 
	Cell* pC;

	for( int k=0; k < number_of_cells ; k++ )
	{ 
		std::vector<double> position = {0,0,0}; 
		position[0] = bounds[0] + UniformRandom()*Xrange; 
		position[1] = bounds[1] + UniformRandom()*Yrange; 
					
		pC = create_cell( *pCD ); 
		pC->assign_position( position );
					
	}

	return; 
}

void save_cells_csv(std::string filename) {

	std::string full_path = PhysiCell_settings.folder + "/" + filename;
    
    // Create and open file
    std::ofstream fileout(full_path);
    
    // Check if file opened successfully
    if (!fileout.is_open()) {
        std::cerr << "Error: Could not create file " << full_path << std::endl;
        return;
    }

	// Write to file
    fileout << "x,y,z,type,volume,cycle entry,custom:GFP,custom:sample" << std::endl;

	for (auto cell : *all_cells) {
		fileout << cell->position[0] << "," << cell->position[1] << "," << cell->position[2] << "," << cell->type_name << std::endl;
	}

	fileout.close();
	std::cout << "Saving cells position at " << full_path << std::endl;

	return;
}

/*************************************/
/*  DIVISION ORIENTATATION FUNCTIONS */
/*************************************/

double cell_neighbor_distance(Cell* pC1, Cell*pC2)
{
	return abs(sqrt(pow((pC1->position[0] - pC2->position[0]), 2) + pow((pC1->position[1] - pC2->position[1]), 2) + pow((pC1->position[2] - pC2->position[2]), 2)));
}

std::vector<Cell*> find_closest_neighbors(const std::vector<Cell*>& cells, Cell* pC) {
		// Using priority queue approach
		auto comp = [pC](Cell* pC1, Cell* pC2) {
			return cell_neighbor_distance(pC1, pC) < cell_neighbor_distance(pC2, pC);
		};
		
		std::priority_queue<Cell*, std::vector<Cell*>, decltype(comp)> pq(comp);
		
		for (Cell* c : cells) {
			pq.push(c);
			if (pq.size() > 3) {
				pq.pop();
			}
		}
		
    std::vector<Cell*> result;
    while (!pq.empty()) {
        result.push_back(pq.top());
        pq.pop();
    }
    std::reverse(result.begin(), result.end());
    return result;
}

//Used instead of UniformOnUnitSphere for cell_division_orientation
std::vector<double> custom_division_orientation( Cell* pC)
{
	if( default_microenvironment_options.simulate_2D == true )
	{		
		std::vector<double> orientation_vec;
		if (pC->type_name == "epi_basal")
		{
			std::vector<Cell*> membrane_neighbors;
			for (Cell* neighbor : pC->state.neighbors)
			{
				if ( neighbor->type_name == "membrane" )
				{ 
					membrane_neighbors.push_back(neighbor); 
				}
			}
			/*std::cout << "There are " << membrane_neighbors.size() << " membrane neighbors." << std::endl;
			std::cout << "Cell in x position: " << pC->position[0] << std::endl;
			std::cout << "In x position: " ;
			for ( int i = 0; i < membrane_neighbors.size(); i++)
			{
			std::cout << membrane_neighbors[i]->position[0] << ", ";
			} 
			std::cout << std::endl;
			*/
			if ( membrane_neighbors.size() > 2)
			{
				//TO TEST FOR 3D, doesn't occur in 3D in our configuration
				std::vector<Cell*> closest_neighbors = find_closest_neighbors(membrane_neighbors, pC);
				std::vector<double> X = {closest_neighbors[0]->position[0], closest_neighbors[1]->position[0], closest_neighbors[2]->position[0]};
				std::vector<double> Y = {closest_neighbors[0]->position[1], closest_neighbors[1]->position[1], closest_neighbors[2]->position[1]};
				std::vector<double> abr2;
				if(X.size() == Y.size())
				{
					abr2 = linreg(X.size(), X, Y);
					std::cout << abr2[0] << ", " << abr2[1] << ", " << abr2[2] << std::endl;
				}
				std::cout << "I am a dum dum this is occuring" << std::endl;
				std::vector<double> result_vec = {-1, abr2[0], 0};
				normalize(&result_vec);
				return result_vec;
			}
			else if ( membrane_neighbors.size() == 2 )
			{
				//This appear to work
				int a = UniformInt();
				if ( a%2 == 0)
				{
					orientation_vec = {membrane_neighbors[0]->position[0] - membrane_neighbors[1]->position[0], membrane_neighbors[0]->position[1] - membrane_neighbors[1]->position[1], membrane_neighbors[0]->position[2] - membrane_neighbors[1]->position[2]};
				}
				else
				{
					orientation_vec = {membrane_neighbors[1]->position[0] - membrane_neighbors[0]->position[0], membrane_neighbors[1]->position[1] - membrane_neighbors[0]->position[1], membrane_neighbors[1]->position[2] - membrane_neighbors[0]->position[2]};
				}
				normalize(&orientation_vec);
				//std::cout << "Orientation vec: " << orientation_vec[0] << "," << orientation_vec[1] << "," << orientation_vec[2] << std::endl;
				//std::cout << "This ran" << std::endl;
				return orientation_vec;
			}
			else if ( membrane_neighbors.size() == 1 )
			{
				std::vector<double> membrane_basal_vector = {membrane_neighbors[0]->position[0] - pC->position[0], membrane_neighbors[0]->position[1] - pC->position[1], membrane_neighbors[0]->position[2] - pC->position[2]};
				// In 2D
				if ( default_microenvironment_options.simulate_2D )
				{ 
					std::vector<double> vect_Z = {0,0,1};
					int a = UniformInt();
					if ( a%2 == 0)
					{
						orientation_vec = cross_product(membrane_basal_vector, vect_Z);
					}
					else
					{
						//not tested yet
						orientation_vec = cross_product(vect_Z, membrane_basal_vector);
					}
					return orientation_vec;
				}
				//In 3D
				else
				{
					//Any vector in the plane perpendicular to membrane_basal_vector will do
					std::vector<double> rand_vec = UniformOnUnitSphere();
					orientation_vec = cross_product(membrane_basal_vector, rand_vec);

					return orientation_vec;
				}
			}
			else
			{
				std::cout << "No membrane neighbor should this epi basal cell really divide ? " << std::endl;
				return UniformOnUnitSphere();
			}
			//Look for for Cells Neighbors of type membrane
			//if 2D, 
				//if more than 2 membrane neighbours
					//draw line going through 3 or 2 closest memebrane cells centers
				//if one membrane neighbour 
					//compute perpendicular axis to cell center neighbour center in xy plane
			
			//if 3D, 
				//if 2 or more neighbour
					//find plane from one segment and 1 point
				//if 1 neighbour
					// Find plane perpendicular to cell membrane cell center

			
		}
		else
		{
			return UniformOnUnitSphere();
		}
	}
	else 
	{
		return UniformOnUnitSphere();
	}
}

/***********************/
/*  LINEAR REGRESSION  */
/***********************/

std::vector<double> linreg(int n, std::vector<double> X, std::vector<double> Y)
{
	double a;
	double b;
	double r;

	double sumx = 0.0;
	double sumx2 = 0.0;
	double sumxy = 0.0;
	double sumy = 0.0;
	double sumy2 = 0.0;

	for (int i=0; i<n; i++)
	{
		sumx += X[i];
		sumx2 += X[i]*X[i];
		sumxy += X[i]*Y[i];
		sumy += Y[i];
		sumy2 += Y[i]*Y[i];
	}
	double denom = (n * sumx2 - sumx*sumx);
	if (denom == 0)
	{
		std::cout << "Error cannot solve the linear regression problem" << std::endl;
		double a = 0.0;
		double b = 0.0;
		double r = 0.0;
	}
	else
	{
		a = (n * sumxy - sumx * sumy)/denom;
		b = (sumy - a * sumx)/n;
		r = (sumxy - sumx * sumy / n)/sqrt((sumx2 - sumx*sumx/n) * (sumy2 - sumy*sumy/n));

		//double a1 = (n * sumxy - sumx * sumy)/denom;
		//double b1 = (sumy * sumx2 - sumxy * sumx)/denom;
		//double r1 = (sumxy - sumx * sumy / n)/sqrt((sumx2 - sumx*sumx/n) * (sumy2 - sumy*sumy/n));
		//std::cout << a1 << ", " << b1 << ", " << r1 << std::endl;
	}
	return {a, b, r*r};
}