#ifndef FLUORESCENCE_CALCULATIONS_H
#define FLUORESCENCE_CALCULATIONS_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "peak_input.h"
#include "constants_input.h"

using namespace std;

double energy_limit = 9.000;
//upper limit for fluorescence integration; limited by CHIANTI database range

int fluorescence_calc(double energy_start, double energy_step, string comet_name, double scaling_factor, double fluo_ratio){

    double u = input_u;

    //inputs intensity spectra and defines interpolation function
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

    gsl_interp_accel *accel_ptr;
    gsl_spline *spline_ptr;
    accel_ptr = gsl_interp_accel_alloc ();
    spline_ptr = gsl_spline_alloc (gsl_interp_cspline, 1000);
    gsl_spline_init (spline_ptr, spectra_energy, spectra_intensity, 1000);


 	//Creates output file and defines number of row in output data based on number of peaks analyzed
	string Inten_name = "../Results/" + comet_name + "/fluorescence_total_output_" + comet_name + ".dat";
	fstream Inten(Inten_name.c_str(), fstream::out | fstream::app);
	Inten << number_of_peaks << endl;
	Inten.close();


    //reads in peak data from peak_input.h
    for( int counter = 1; counter <= number_of_peaks; counter++ ){

        string peak_name = input_peak_name[counter - 1];
        double peak_energy = input_peak_energy[counter - 1];
        double peak_efficiency = input_peak_efficiency[counter - 1];
        double element_mass = input_element_mass[counter - 1];

        if (peak_energy <= energy_start || peak_energy >= energy_limit){
            cout << "Warning: Fluorescence peak " << peak_name << " was not included in the model as it lies outside the designated energy range." << endl;
            cout << "Please adjust either the energy range or seleced fluorescence peaks in future runs. " << endl;
        }

        else{
            //multiplies ratio with efficiency rate
            double peak_ratio = fluo_ratio * peak_efficiency;

            //inputs the cross-sections from a .dat file
            string input_name = "../Inputs/" + peak_name + "_cross.dat";
            fstream input_file(input_name.c_str(), fstream::in);

            int row;
            input_file >> row;

            double input[row][2];
            for( int i=0; i<row; i++){
                for( int j=0; j<2; j++){
                    input_file >> input[i][j]; } }
            input_file.close();

            double energy[row];
            double sigma[row];
            double ratio[row];
            for (int l=0; l<row; l++){
                energy[l] = (input[l][0]); //in keV

                sigma[l] = input[l][1]/(1e4) * element_mass * u;
                //converts total cross-section from cm^2/g to m^2/g
                //and multiplies by atomic mass

                ratio[l] = peak_ratio*sigma[l];
            }

            gsl_interp_accel *accel_ptr_ele;
            gsl_spline *spline_ptr_ele;
            accel_ptr_ele = gsl_interp_accel_alloc ();
            spline_ptr_ele = gsl_spline_alloc (gsl_interp_cspline, row);
            gsl_spline_init (spline_ptr_ele, energy, ratio, row);

            //calculates total intensity from fluorescence for a single peak
            double total_integral = 0;
            double energy_temp = peak_energy;
            while (energy_temp <= energy_limit){

                double intensity_temp = gsl_spline_eval (spline_ptr, energy_temp, accel_ptr);
                double ratio_temp = gsl_spline_eval (spline_ptr_ele, energy_temp, accel_ptr_ele);

                total_integral = total_integral + intensity_temp * energy_step * ratio_temp;

                energy_temp = energy_temp + energy_step;
            }

            //outputs fluorescence components to file
            fstream Inten(Inten_name.c_str(), fstream::out | fstream::app);
            Inten << scientific << peak_energy << " " << total_integral << endl;
            Inten.close();

            gsl_spline_free (spline_ptr_ele);
            gsl_interp_accel_free (accel_ptr_ele);
            }
    }

    gsl_spline_free (spline_ptr);
    gsl_interp_accel_free (accel_ptr);

    return 0;
}

#endif
