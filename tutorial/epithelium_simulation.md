# Simulating epithelium 101

**Navigate the project of epithelium modelling currently named “test_start_and_stop”:**

1. Download the version of PhysiCell I uploaded in GitLab named epi_model
2. Then go to this PhysiCell repo and run the following command:

```bash
make load PROJ=test_start_and_stop
make #generate an exe file named heterogeneity
./heterogeneity
```

This will build the PhysiCell project for epithelium modelling and then run a simulation with the default parameters for this project. The output files (.svg, .txt, .mat, …) are saved by default on the ./PhysiCell/output/ folder.

To modifiy the parameters of the simulation, you will need to modify the file PhysiCell_settings.xml in the /PhysiCell/config/ folder. You also have the possibility to generate an other .xml file (you can use PhysiCell Studio to do so) and then to run your program using this parameter file run the command:

```bash
./heterogeneity ./path_to_your_file_from_PhysiCell_folder/New_settings_file.xml
```

You can also modify the cells.csv file or create a new one (then specify the new name of file in the .xml), this file indicates the position and type of cells at t=0. It is also possible to generate them by runnning ./heterogeneity, see next section for more information.

It is as well possible to modify or create a new cell_rules.csv (don’t forget to then specify the name change in your .xml). This way you can alter how cells interact with their environment, and with each other. 

**Success you can now run the model as you wish !**

1. For more advanced modifications of the project: The majority of what defines the model is specified in the files custom.cpp and custom.h of the project. If you modify those files in the test_start_and_stop folder you will need to execute again the commands in 2. for the modifications to be taken into account. Careful if you modify the files in the ./PhysiCell/custom_modules folder, you will just need to run `make`  again. But these will be erased if you run the commands in 4.
2. How to start fresh ? (clean folder and go back to default): You can run the following commands to erase files in output folder, erase the .o files and go back to PhysiCell defaults:

```bash
make data-cleanup #erases all files in ./output
make clean #erases all .O files in ./
make reset #erases main.cpp, files in ./config, and in ./custom_modules and replaces with default files
```

By running these 3 commands you go back to the configuration of PhysiCell you downloaded! If you made modifications in ./user_projects/test_start_and_stop or in an otehr uer project they will not be affected. However, other modifications will most likely be erased!! So make sure you saved everything you needed before running these.

**For more information about how to start using PhysiCell the file [Quickstart.md](http://Quickstart.md) in the folder documentation-deprecated can still be useful!**

**Generate files for initial cells position:**

To generate one or multiple initial epithelium file (cells.csv files), you can when you have generated your `heterogeneity` file run the command:

```bash
./heterogeneity ./path_to_file/New_settings_file.xml 1 default
```

The first argument  is the path to the config file you want to use, for instance if you want the default write: `./config/PhysiCell_settings.xml` . 

The second argument is the number of file you want to generate. They will be saved in the `./output` folder and be named cells_0.csv, cells_1.csv and so on. 

The third argument for now is a fraud, whatever you type the same kind of epithelium will be generated, BUT you can add your own function in custom.cpp and then create anything you want !

**Use the “start and stop” to interrupt a simulation when a condition is met and start a simulation from saved files:**

In custom.cpp you can write your own auto stop function and declare it in custom.h. How to proceed in detail to add a new condition for autostop: 

In custom.cpp: 

1. Create an unordered map that associates the user_parameters related to start and stop and a boolean value that indicates the value if a bool, or if the value is defined if the parameter is a string. They are initialized to false and then the value are set from the .xml, therefore if a value is not indicated in .xml, this prevents a segmentation fault during the simulation. The function evaluate_start_stop_parameters checks that the parameters are defined in the .xml and attributes the value. It checks as well that the path given exist. You will have to add the parameter relative to your new auto_stop function in auto_stop_param (not actually necessary but encouraged for robust parsing).
    
    ```cpp
    /****************************************/
    /* START AND STOP FUNCTIONS DEFINITIONS */
    /****************************************/
    using namespace std;
    vector<double> vector_alives;
    std::unordered_map<std::string, bool> auto_stop_param = {{"start_stop", false}, {"read_init", false},{"saving_folder", false}, {"init_cells_filename", false}, {"auto_stop", false}, {"auto_stop_alive", false}, {"auto_stop_epi_stable", false}, {"auto_stop_epi_size", false}};
    
    int evaluate_start_stop_parameters() {    
    	int result = 1;    //auto stop parameters should be boolean or string in the user_parameters section of your .xml config files    
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
    ```
    
    **You will need to add evaluate_start_stop_parameters at the end of your setup_microenvironment function.**
    
2. Write your auto_stop function it can take any arguments as input but it has to return a bool, false if the condition is not met and true otherwise, for instance the following function aims to stop the simulation once  the number of alive cells is stable:

```cpp
bool auto_stop_alive(int alive_cells) {    //concatenate the number of alive cells to the vector    
	vector_alives.push_back(alive_cells);    
	std::cout << "Steps: " << vector_alives.size() << std::endl;
  bool condition = false;    // check the number of elements inside the vector to decide if process it and compute the derivative    
  if (vector_alives.size() >= 8) {        
	  std::vector<double> derivative; // compute the derivative only for the last three steps        
    for (size_t i = vector_alives.size() - 4; i < vector_alives.size(); ++i) {            
	    double slope = vector_alives[i] - vector_alives[i - 1];            
	    derivative.push_back(slope);        
		}
    condition = true;        
    for (double slope : derivative) {            // if the slope is less than or equal to zero, set condition to false            
	    if (slope > 100) {               
		    condition = false;                
		    break;            
		  }        
	  }
  } 
  else { condition = false; }
  
  bool stop;
  if (condition) { stop = true; } 
  else { stop = false; }    
  return stop;
 }
```

Then you are going to add a condition to compute this function at the appropriate time in your main.cpp:

```cpp
// INSERT HERE YOUR AUTO STOP FUNCTION                    
//These conditions are evaluated only at full_save times                    
if(parameters.bools("auto_stop")){                        
	int alive = total_live_cell_count();                        
	std::vector<double> vector_epi_size;
                        
  //Computation necessary for both auto_stop_epi_size and auto_stop_epi_stable                        
  if(parameters.bools("auto_stop_epi_size") || parameters.bools("auto_stop_epi_stable")){                            
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
  if(parameters.bools("auto_stop_epi_size")){                                
	  stop = auto_stop_epi_size(vector_epi_size, parameters.doubles("epi_max_size"));                        
	}
	
  //auto stop condition stable epi size                        
  if(parameters.bools("auto_stop_epi_stable")){                            
	  //Only check after a certain time                            
	  if (PhysiCell_globals.current_time > PhysiCell_settings.full_save_interval * 20)                            
	  {                                
		  stop = auto_stop_epi_stable(vector_epi_size, deque_epi_average_size, 10, 3);                            
		}                        
	}                    
}
```

In your .xml: As said previously, you will need to add a parameter in your .xml file, this will allow to enable/disable this condition from you .xml file, you can write in the .xml as such:

```xml
<user_parameters>
	<auto_stop_epi_size type="bool" units="">false</auto_stop_epi_size>
</user_parameters>
```

Then you can access this parameter in your main.cpp using :

```cpp
parameters.bools("auto_stop_epi_size")
```

Or with :

```cpp
auto_stop_param["auto_stop_epi_size"]
```

If you have added to the unordered map auto_stop_param `{”auto_stop_epi_size”, false}` .

You will most likely need to create new objects to provide inputs for your auto_stop function, just be careful where you put them and when they have to be computed!