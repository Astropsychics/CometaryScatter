#ifndef ATOMIC_CROSS_SECTION_INPUT_H
#define ATOMIC_CROSS_SECTION_INPUT_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "constants_input.h"

using namespace std;

int atomic_cross(double energy_start, double energy_end, double energy_step, string comet_name, double scaling_factor, double ratio, double Brightness){

    //input from constant_input.h
	double u = input_u;

	//inputs intensity spectra and generates interpolation function
	fstream intensity_file("../Inputs/Spectra/Total_Intensity_CHIANTI.dat", fstream::in);
    double intensity_input[1000][2];
    for( int i=0; i<1000; i++ ){
        for( int j=0; j<2; j++ ){
            intensity_file >> intensity_input[i][j]; } }
    intensity_file.close();

    double spectra_energy[1000];
    double spectra_intensity[1000];
    for (int l=0; l<1000; l++){

    spectra_energy[l] = intensity_input[l][0];
    spectra_intensity[l] = scaling_factor * intensity_input[l][1]; }

    gsl_interp_accel *accel_ptr_inten;
    gsl_spline *spline_ptr_inten;
    accel_ptr_inten = gsl_interp_accel_alloc ();
    spline_ptr_inten = gsl_spline_alloc (gsl_interp_cspline, 1000);
    gsl_spline_init (spline_ptr_inten, spectra_energy, spectra_intensity, 1000);


	//input parameters for elements
	string element[4] = {"H","C","N","O"};
	double mass[4] = {1.00794, 12.0107, 14.0071, 15.9994};

    //inputs cross sections
	for (int x = 0; x < 4; x++){
		string input_name = "../Inputs/Atomic_Cross_Sections/" + element[x] + "_sigma.txt";
		ifstream input_file(input_name.c_str());

  		int row;
        input_file >> row;

		double input[row][2];
  			for( int i=0; i<row; i++){
        		for( int j=0; j<2; j++){
        	    	input_file >> input[i][j]; } }
        input_file.close();

		double energy[row];
		double sigma[row];
 		for (int l=0; l<row; l++){

			energy[l] = (input[l][0]);  // in keV

  			sigma[l] = input[l][1]/(1e4) * mass[x] * u; }
			//converts total cross-section from cm^2/g to m^2/g and multiplies by atomic mass
			//to convert it to a differential cross section in units of m^2

   		gsl_interp_accel *accel_ptr;
   		gsl_spline *spline_ptr;
  		accel_ptr = gsl_interp_accel_alloc ();
  		spline_ptr = gsl_spline_alloc (gsl_interp_cspline, row);
  		gsl_spline_init (spline_ptr, energy, sigma, row);



		//multiplies ratio with cross section and outputs to file
		string final_name = "../Results/" + comet_name + "/" + element[x] + "_output_" + comet_name + ".dat";
		ofstream output(final_name.c_str());

 		double energy_temp = energy_start;
		while (energy_temp <= energy_end){

			double solar_intensity_temp = gsl_spline_eval (spline_ptr_inten, energy_temp, accel_ptr_inten);
			double cross_temp = gsl_spline_eval (spline_ptr, energy_temp, accel_ptr);
			double Intensity = ratio * cross_temp * solar_intensity_temp;

			output << scientific << energy_temp << " " << cross_temp
			<< " " << Intensity << " " << Brightness*Intensity << endl;
			energy_temp = energy_temp + energy_step;
        }

		output.close();
		gsl_spline_free (spline_ptr);
	   	gsl_interp_accel_free (accel_ptr);

	}

	gsl_spline_free (spline_ptr_inten);
   	gsl_interp_accel_free (accel_ptr_inten);

	return 0;
}

#endif
