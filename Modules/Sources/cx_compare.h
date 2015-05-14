#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "../comet_input.h"

using namespace std;

int cx_compare(double energy_start, double energy_end, double energy_step, int comet_number){

    //inputs comet name from comet_input.h
    string comet_name = input_comet_name[comet_number - 1];

    //calculates numbers of steps and creates appropriately-sized array
    double energy = energy_start;
    int energy_row = (energy_end - energy_start)/energy_step + 1;

    //input observation file
    string observation_name = "../Inputs/Observations/" + comet_name + "_intensity_Chandra.dat";
    ifstream observation_file(observation_name.c_str());

    int obs_row;
    observation_file >> obs_row;

    double input_obs[obs_row][2];
    double obs_energy[obs_row];
    double obs_intensity[obs_row];
    for( int i=0; i<obs_row; i++ ){
        for( int j=0; j<2; j++ ){
            observation_file >> input_obs[i][j]; }
        obs_energy[i] = input_obs[i][0];
        obs_intensity[i] = input_obs[i][1];
    }
	observation_file.close();

    gsl_interp_accel *accel_ptr_obs;
    gsl_spline *spline_ptr_obs;
    accel_ptr_obs = gsl_interp_accel_alloc ();
    spline_ptr_obs = gsl_spline_alloc (gsl_interp_cspline, obs_row);
    gsl_spline_init (spline_ptr_obs, obs_energy, obs_intensity, obs_row);

    //input model file
    fstream model_file("../Results/CX_Chandra_spectrum.dat", fstream::in);

    double input_model[energy_row][2];
    double model_energy[energy_row];
    double model_intensity[energy_row];
    for( int i=0; i<energy_row; i++ ){
        for( int j=0; j<2; j++ ){
            model_file >> input_model[i][j]; }
        model_energy[i] = input_model[i][0];
        model_intensity[i] = input_model[i][1];
    }
    model_file.close();

    gsl_interp_accel *accel_ptr_model;
    gsl_spline *spline_ptr_model;
    accel_ptr_model = gsl_interp_accel_alloc ();
    spline_ptr_model = gsl_spline_alloc (gsl_interp_cspline, energy_row);
    gsl_spline_init (spline_ptr_model, model_energy, model_intensity, energy_row);


    //calculates the average ratio between the observational data and the model to find an
    //optimal scaling factor
    double scaling_factor = 0;
    double counter = 0;

    for ( float energy = 0.400 ; energy <= 0.750; energy += 0.010 ){
        double obs_temp = gsl_spline_eval (spline_ptr_obs, energy, accel_ptr_obs);
        double model_temp = gsl_spline_eval (spline_ptr_model, energy, accel_ptr_model);

        scaling_factor += obs_temp / model_temp;
        counter++;
    }

    double average_scaling_factor = scaling_factor/counter;
    
    //outputs re-scaled Chandra spectrum
    string output_name_Chandra = "../Results/" + comet_name + "/Chandra_CX_spectrum_"+ comet_name + ".dat";
    ofstream output_file_Chandra(output_name_Chandra.c_str());
    for( int i=0; i<energy_row; i++ ){
        output_file_Chandra << scientific << model_energy[i] << " " << model_intensity[i]*average_scaling_factor << endl;
    }
    output_file_Chandra.close();

    //inputs intensity spectrum
    fstream spectrum_file("../Results/CX_spectrum.dat", fstream::in);
    double input_spectrum[energy_row][2];
    double spectrum_energy[energy_row];
    double spectrum_intensity[energy_row];
    for( int i=0; i<energy_row; i++ ){
        for( int j=0; j<2; j++ ){
            spectrum_file >> input_spectrum[i][j]; }
        spectrum_energy[i] = input_spectrum[i][0];
        spectrum_intensity[i] = input_spectrum[i][1];
    }
    spectrum_file.close();

    //outputs re-scaled intensity spectrum
    string output_name_spectrum = "../Results/" + comet_name + "/CX_spectrum_"+ comet_name + ".dat";
    ofstream output_file_spectrum(output_name_spectrum.c_str());
    for( int i=0; i<energy_row; i++ ){
        output_file_spectrum << scientific << spectrum_energy[i] << " " << spectrum_intensity[i]*average_scaling_factor << endl;
    }
    output_file_spectrum.close();


    gsl_spline_free (spline_ptr_obs);
    gsl_interp_accel_free (accel_ptr_obs);

    gsl_spline_free (spline_ptr_model);
    gsl_interp_accel_free (accel_ptr_model);

    return 0;
}
