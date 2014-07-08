#ifndef DUST_MIXING_RATIO_H
#define DUST_MIXING_RATIO_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>

using namespace std;

int dust_mixing_ratio(double energy_start, double energy_end, double energy_step, string comet_name){

	//input parameters for elements
	string element[3] = {"C_dust","H2O","Si_dust"};
	double mixing_ratio[3] = {0.1, 0.85, 0.05};
	
	int row = (energy_end - energy_start)/energy_step + 1;
	vector< vector<double> > total_intensity(row, std::vector<double>(3, 0));

	for (int x=0; x < 3; x++) {
	
		string input_name = "../Results/" + comet_name + "/" + element[x] + "_output_" + comet_name + ".dat";
		ifstream input_file(input_name.c_str());
    	
		vector< vector<double> > input(row, std::vector<double>(3, 0));
    		for( int i=0; i<row; i++){
        		for( int j=0; j<3; j++){
        	    		input_file >> input[i][j]; } }
    		input_file.close();


		for( int k=0; k<row; k++){
			total_intensity[k][0] = input[k][0];
			total_intensity[k][1] = total_intensity[k][1] + mixing_ratio[x] * input[k][1];
			total_intensity[k][2] = total_intensity[k][2] + mixing_ratio[x] * input[k][2]; 
		}	
	}

	string final_name = "../Results/" + comet_name + "/dust_total_output_" + comet_name + ".dat";
	ofstream output(final_name.c_str());
 	
	for ( int l=0; l<row; l++ ){
		output << scientific << total_intensity[l][0]  << " " << total_intensity[l][1]
		<< " " << total_intensity[l][2] << endl; } 
	output.close();

	return 0;
}

#endif
