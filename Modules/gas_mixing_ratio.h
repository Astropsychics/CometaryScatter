#ifndef GAS_MIXING_RATIO_H
#define GAS_MIXING_RATIO_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>

#include "comet_input.h"

using namespace std;

int gas_mixing_ratio(double energy_start, double energy_end, double energy_step, string comet_name){

    //inputs number of different gases included from comet_input.h
    int gas_types = input_gas_types;
	
    //determines number of energy steps to perform and generates intensity array
	int row = (energy_end - energy_start)/energy_step + 1;
	vector< vector<double> > total_intensity(row, std::vector<double>(3, 0));
	
	for (int x = 0; x < gas_types; x++){
        
        string element = input_gas_element[x];
        double mixing_ratio = input_gas_mixing_ratio[x];
        
		//reads in gas spectrum files
        string input_name = "../Results/" + comet_name + "/" + element + "_output_" + comet_name + ".dat";
		ifstream input_file(input_name.c_str());
    	
		vector< vector<double> > input(row, std::vector<double>(4, 0));
    		for( int i=0; i<row; i++){
        		for( int j=0; j<4; j++){
        	    	input_file >> input[i][j]; } }
    		input_file.close();

        //applies mixing ratios
		for( int k=0; k<row; k++){
			total_intensity[k][0] = input[k][0];
			total_intensity[k][1] = total_intensity[k][1] + mixing_ratio * input[k][2];
			total_intensity[k][2] = total_intensity[k][2] + mixing_ratio * input[k][3]; }
	}
    
    //writes to output file
	string final_name = "../Results/" + comet_name + "/gas_spectrum_" + comet_name + ".dat";
	ofstream output(final_name.c_str());
		
 	for ( int l=0; l<row; l++ ){
		output << scientific << total_intensity[l][0]  << " " << total_intensity[l][1] 
			<< " " << total_intensity[l][2] << endl; } 
	output.close();

	return 0;
}

#endif
