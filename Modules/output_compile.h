#ifndef OUTPUT_COMPILE_H
#define OUTPUT_COMPILE_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <sstream>
#include <vector>

using namespace std;

//Comet info
double grain_size_upper = 99.0; 

int output_compile(double energy_start, double energy_end, double energy_step, string comet_name, int x ){
	
	string element[3] = {"C","H2O","Si"};
	string output_element[3] = {"C_dust","H2O","Si_dust"}; 

	int row = (energy_end - energy_start)/energy_step + 1;

	double energy_array[row]; 
	vector<double> compile_intensity_array(row,0);
	vector<double> compile_brightness_array(row,0);
    
    	for( int a=11; a <= grain_size_upper; a++ ){

		string element_name = element[x];
		//additional declaration of element name required due to 
		//ostringstream conversion as it converts grain_size to a 
		//const char*, which cannot be added to a const char [x]
		//when defining the input file name

	
		double grain_size = a/10.0;
		string grain_string = static_cast<ostringstream*>( &(ostringstream() << grain_size) )->str();

        	string input_name = "../Results/" + element_name + "_Results/" + comet_name + "/" 
		+ element_name + "_output_" + grain_string + "nm_" + comet_name + ".dat";
        	fstream input_file(input_name.c_str(), fstream::in);

		vector< vector<double> > ratio_input(row, std::vector<double>(3, 0));
    		for( int i=0; i<row; i++ ){
       	    		for( int j=0; j<3; j++ ){
                		input_file >> ratio_input[i][j]; } }
        	input_file.close();

    		for( int k=0; k<row; k++ ){
	    		energy_array[k] = ratio_input[k][0];

           		compile_intensity_array[k] = compile_intensity_array[k] + ratio_input[k][1];
	    		compile_brightness_array[k] = compile_brightness_array[k] + ratio_input[k][2];
          
    		}
    	} 

	string final_name = "../Results/" + comet_name + "/" + output_element[x] + "_output_" + comet_name + ".dat";
    	fstream output(final_name.c_str(), fstream::out);
    
    	for( int l=0; l<row; l++ ){ 
		output << scientific << energy_array[l] << " " << compile_intensity_array[l] 
		<< " " << compile_brightness_array[l] << endl;
	}
    	output.close();

 	return 0;
}

#endif

