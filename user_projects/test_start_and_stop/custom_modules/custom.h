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
#include <sstream>
#include <vector>
#include <string>
#include <deque>

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
#include "../addons/start_and_stop/start_and_stop.h"

using namespace BioFVM; 
using namespace PhysiCell;

// setup functions to help us along 

void create_cell_types( void );
void setup_tissue( void );
void setup_tissue_input(void); //Necessary to restrat simulation from input

// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 

// custom pathology coloring function 

std::vector<std::string> my_coloring_function( Cell* );

/**********************/ 
/* CUSTOM FUNCTIONS   */ 
/**********************/ 


//functions for initialization

// helper function to read init files
std::vector<std::vector<double>>  read_cells_positions(std::string filename, char delimiter, bool header);

// helper function to create a sphere of cells of a given radius
std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius);

// helper function to create a disc of cells of a given radius
std::vector<std::vector<double>> create_cell_disc_positions(double cell_radius, double disc_radius);

// helper function that calculates phere volume
inline float sphere_volume_from_radius(float radius) {return 4/3 * PhysiCell_constants::pi * std::pow(radius, 3);}

// helper function to inject density surrounding a spheroid
void inject_density_sphere(int density_index, double concentration, double membrane_lenght);

// helper function to remove a density
void remove_density( int density_index );

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt );
void custom_function( Cell* pCell, Phenotype& phenotype , double dt );

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt ); 

std::vector<std::string> heterogeneity_coloring_function( Cell* );

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt ); 
double total_live_cell_count();

// count the number of total dead cells at current time step
double total_dead_cell_count();

// function to auto stop the simulation
bool auto_stop_resistance(int alive_cells, int resistant_cells);    
bool auto_stop_epi_size(std::vector<double> vector_epi_pos, double epi_max_size);
bool auto_stop_epi_stable(std::vector<double> vector_epi_pos, std::deque<double> &deque_epi_average_size, double steps, double tolerance);
bool auto_stop_alive(int alive_cells);
int save_resistant_cells(std::ofstream& file_resistant);
bool auto_stop();

//Functions to create an pre-epithelium
void create_pre_epithelium(int argc, char* argv[]);
void create_epithelium_csv(int nbr, std::string func);
void position_epithelium_cells(void);
void random_fill_rectangle (BioFVM::gradient bounds, PhysiCell::Cell_Definition *pCD, double confluence);
void save_cells_csv(std::string filename);

//Functions for divisions orientation
double cell_neighbor_distance(Cell* pC1, Cell*pC2);
std::vector<Cell*> find_closest_neighbors(const std::vector<Cell*>& cells, Cell* pC); // Not tested yet

std::vector<double> custom_division_orientation(Cell* pC);
