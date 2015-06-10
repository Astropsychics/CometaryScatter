#ifndef DUST_LARGE_GRAIN_H
#define DUST_LARGE_GRAIN_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>

#include "./dust_input.h"

using namespace std;

int dust_large_grain(double energy_start, double energy_end, double energy_step, string comet_name){

    //inputs dust model variables from dust_input.h
    double alpha = input_alpha; //Dust density dependence on grain radius
    double gamma = input_gamma; //Dust density dependence on cross section
    double a_max = input_a_max; //Maximum dust radius

    //determines numbers of energy steps to perform and generates intensity array
    int row = (energy_end - energy_start)/energy_step + 1;
    vector< vector<double> > total_intensity(row, std::vector<double>(3, 0));

    //Calculates normalization constant.
    //Step size of integration is 1e-10
    double norm_const = 0;
    for ( double a = 1e-9; a <= 1e-8; a += 1e-10){
        norm_const += pow(a, gamma - alpha) * 1e-10;
    }

    //Calculates additional contribution from larger grains.
    //Step size of integration is 1e-8
    double large_grain_int = 0;
    for ( double a = 1e-8; a <= a_max; a += 1e-8){
        large_grain_int += pow(a, gamma - alpha) * 1e-8;
    }

    //inputs dust spectrum file
	string input_name = "../Results/" + comet_name + "/dust_spectrum_small_grain_" + comet_name + ".dat";
	ifstream input_file(input_name.c_str());

	vector< vector<double> > input(row, std::vector<double>(3, 0));
    	for( int i=0; i<row; i++){
        	for( int j=0; j<3; j++){
         		input_file >> input[i][j]; } }
    	input_file.close();



	for( int k=0; k<row; k++){
		total_intensity[k][0] = input[k][0];
		total_intensity[k][1] = input[k][1] * (1 + large_grain_int / norm_const);
		total_intensity[k][2] = input[k][2] * (1 + large_grain_int / norm_const);
	}


    //Write output to file
	string final_name = "../Results/" + comet_name + "/dust_spectrum_" + comet_name + ".dat";
	ofstream output(final_name.c_str());

	for ( int l=0; l<row; l++ ){
		output << scientific << total_intensity[l][0]  << " " << total_intensity[l][1]
		<< " " << total_intensity[l][2] << endl; }
	output.close();

	return 0;
}

#endif
