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
# Copyright (c) 2015-2022, Paul Macklin and the PhysiCell Project             #
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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>
#include <sys/stat.h>
#include <functional>
#include <algorithm>
#include <deque>

#include "./core/PhysiCell.h"
#include "./core/PhysiCell_cell.h"
#include "./core/PhysiCell_utilities.h"
#include "./modules/PhysiCell_settings.h"
#include "./modules/PhysiCell_standard_modules.h" 

// put custom code modules here! 

#include "./custom_modules/custom.h" 
	
using namespace BioFVM;
using namespace PhysiCell;

int main( int argc, char* argv[] )
{
	//Declaring the .xml file path
	std::string xml_path_str;
	
	bool XML_status = false; 
	char copy_command[1024]; 
	
	if( argc > 3 )
	{
		create_pre_epithelium(argc, argv);

		sprintf( copy_command , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() );
		system( copy_command );

		char filename[1024];
		sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 

		save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , 0.0 ); 
		// save a quick SVG cross section through z = 0, after setting its 
		// length bar to 200 microns 

		PhysiCell_SVG_options.length_bar = 200; 

		// for simplicity, set a pathology coloring function 
		
		std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function; 
		sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
		SVG_plot( filename , microenvironment, 0.0 , 0.0 , cell_coloring_function, -50.0);
		
		sprintf( filename , "%s/legend.svg" , PhysiCell_settings.folder.c_str() ); 
		create_plot_legend( filename , cell_coloring_function ); 


		return 0;
		
	}

	// Keeping track of time when starting stopping and reloading the simulation
	clock_t T_save_start, T_save_stop, T_reload_start, T_reload_stop, T_total_start, T_total_stop, T_main_start, T_main_stop;
	T_total_start = clock();
	T_reload_start = clock();

	std::ofstream file_times("output/interesting_times.txt", std::ios::app);

	// load and parse settings file(s)
	if( argc > 1 )
	{
		xml_path_str = argv[1];
		std::cout << argv[1] << std::endl;
		XML_status = load_PhysiCell_config_file( argv[1] ); 
		sprintf( copy_command , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() ); 
	}
	else
	{
		xml_path_str = "./config/PhysiCell_settings.xml"; 
		XML_status = load_PhysiCell_config_file( xml_path_str );
		sprintf( copy_command , "cp ./config/PhysiCell_settings.xml %s" , PhysiCell_settings.folder.c_str() ); 
	}
	if( !XML_status )

		{ exit(-1); }

	// copy config file to output directry 
	system( copy_command ); 
	
	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);

	// PNRG setup (Following line added in main.cpp of example project in start and stop folder)
	//SeedRandom(); // or specify a seed here
	
	// time setup 
	std::string time_units = "min"; 

	/* Microenvironment setup */ 
	setup_microenvironment(); // modify this in the custom code 
	
	// Getting the value of start_stop from user_parameters in the .xml setting file
	std::cout << "Start and stop value: " << auto_stop_param["start_stop"] << std::endl;

	//User Parameters
		//Additional parameters in the setting .xml file that will be used as condition to stop the simulation
	//double stop_time = parameters.doubles("stop_time");
	//double epi_total_thickness = parameters.doubles("max_thickness");
	
	/* PhysiCell setup */ 
 		// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );

	/* Users typically start modifying here. START USERMODS */ 
	create_cell_types();
	if( auto_stop_param["start_stop"] ) {
		auto_stop_param["read_init"] = true;
		// reset cells as they were in the previous simulation
		setup_tissue();
		reset_cell( cell_container->last_cell_cycle_time, parameters.strings("saving_folder"), xml_path_str );
		//exit(-1);

		reset_global_parameters( cell_container, parameters.strings("saving_folder") );
		reset_microenv( parameters.strings("saving_folder") );

	} else{
		setup_tissue(); //death model index = 1 == necrotic...= 0 == apoptotic.
	}

	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = my_coloring_function; 
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, -50.0 );
	
	sprintf( filename , "%s/legend.svg" , PhysiCell_settings.folder.c_str() ); 
	create_plot_legend( filename , cell_coloring_function ); 
	
	display_citations(); 

	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	T_reload_stop = clock();
	
	std::ofstream report_file;
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}

	//put here reset randomness
	if( auto_stop_param["start_stop"] )
	{
		reset_randomness( parameters.strings("saving_folder") );
	}
	if (auto_stop_param["auto_stop"])
	{
		mkdir(parameters.strings("saving_folder").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
	//define auto stop variable
	bool stop = false;
	std::deque<double> deque_epi_average_size;
	bool first_save = false;

	T_main_start = clock();
	
	// main loop 
	
	try 
	{		
		while( PhysiCell_globals.current_time < PhysiCell_settings.max_time + 0.1*diffusion_dt && stop!=true)
		{
			
			//Save at time of first mechanical step
			if( PhysiCell_globals.current_time > phenotype_dt && first_save == false)
			{
				sprintf( filename , "%s/first_save", PhysiCell_settings.folder.c_str());
				save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time );
				first_save = true;
			}
			// save data if it's time. 
			if( PhysiCell_globals.current_time > PhysiCell_globals.next_full_save_time - 0.5 * diffusion_dt )
			{
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
					// INSERT HERE YOUR AUTO STOP FUNCTION
					//These conditions are evaluated only at full_save times
					if(auto_stop_param["auto_stop"]) {
						int alive = total_live_cell_count();
						std::vector<double> vector_epi_size;
						//Computation necessary for both auto_stop_epi_size and auto_stop_epi_stable
						if(auto_stop_param["auto_stop_epi_size"] || auto_stop_param["auto_stop_epi_stable"]){
							size_t n = std::min((*all_cells).size(), size_t(50));
							std::vector<double> y_positions;
							y_positions.reserve((*all_cells).size());
							//auto stop condition (alive)
							for (Cell *cell : *all_cells){
								y_positions.push_back(cell->position[1]);
							}
							// Partially sort to get top 100 largest elements
							std::partial_sort(y_positions.begin(), y_positions.begin() + n, y_positions.end(), std::greater<double>());
							vector_epi_size = std::vector<double>(y_positions.begin(), y_positions.begin() + n);
							
						}
						//auto stop when the epithelium reaches a given size
						if(auto_stop_param["auto_stop_epi_size"]){
							std::cout << "Vector epi size: " << vector_epi_size.size() << std::endl;
							stop = auto_stop_epi_size(vector_epi_size, parameters.doubles("epi_max_size"));
						}
						//auto stop condition stable epi size
						if(auto_stop_param["auto_stop_epi_stable"]){
							//Only check after a certain time
							if (PhysiCell_globals.current_time > PhysiCell_settings.full_save_interval * 20)
							{
								stop = auto_stop_epi_stable(vector_epi_size, deque_epi_average_size, 10, 3);
							}
						}
					}
				}
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			if( PhysiCell_globals.current_time > PhysiCell_globals.next_SVG_save_time - 0.5 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, -50.0 );
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			
			/*
			  Custom add-ons could potentially go here. 
			*/
			
			PhysiCell_globals.current_time += diffusion_dt;
		}
		
		if( PhysiCell_settings.enable_legacy_saves == true )
		{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	T_main_stop = clock();
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, -50.0 );
	
	// Save all the files needed for Start & Stop at the right point.
	
	if (auto_stop_param["auto_stop"]) 
	{
		T_save_start = clock();
		save_cell_microenv_data(cell_container, parameters.strings("saving_folder"));
		std::cout << "cells data saved successfully" << std::endl;
		T_save_stop = clock();
	}

	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	T_total_stop = clock();
	double T_save, T_reload, T_total, T_main;
	T_save = (double)(T_save_stop - T_save_start)/CLOCKS_PER_SEC;
	T_reload = (double)(T_reload_stop - T_reload_start)/CLOCKS_PER_SEC;
	T_total = (double)(T_total_stop - T_total_start)/CLOCKS_PER_SEC;
	T_main = (double)(T_main_stop - T_main_start)/CLOCKS_PER_SEC;
	file_times << T_save << " " << T_reload << " " << T_total <<" " << T_main <<  std::endl;
	file_times.close();
	return 0; 
}
