#include "mpi.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "../comet_input.h"
#include "../constants_input.h"
#include "../task_division.h"
#include "../output_compile.h"
#include "../dust_mixing_ratio.h"

using namespace std;

//Dust model variables
double alpha = 2.5; 		       //Dust density dependence on grain radius
double beta = 2.0; 		       //Dust density dependence on radius from cometary center
double a_low = 1e-9;		       //Lower bound of grain size integration
double delta_a = 0.1e-9;	       //Grain radius step size


int dust_core(double energy_start, double energy_end, double energy_step, int comet_number){

	string element[3] = {"C","H2O","Si"};
		
	int rank;  	
	MPI_Comm_rank( MPI_COMM_WORLD , &rank );

	double Ro = input_Ro;
	double Q = input_Q[comet_number - 1]; 
	double vel = input_vel[comet_number - 1];
	double r = input_r[comet_number - 1];
	double rd = input_rd[comet_number - 1];
	double rg = input_rg[comet_number - 1];
	double scaling_factor = input_scaling_factor[comet_number - 1];
	string comet_name = input_comet_name[comet_number - 1];

	/////////////////////////////////
    	double Np = (Q/4.76e27)*1e27;      //Number of particles, as determined from Fink & Rubin 2012, 
	/////////////////////////////////

	//inputs intensity spectra/////////////////////////////////////////////////////
	fstream intensity_file("../Inputs/Spectra/Total_Intensity_CHIANTI.dat", fstream::in); 
    	double intensity_input[1000][2] = {0};
    	for( int i=0; i<1000; i++ ){
       		for( int j=0; j<2; j++ ){
            		intensity_file >> intensity_input[i][j]; } }
    	intensity_file.close();

    	double spectra_energy[1000] = {0};
    	double spectra_intensity[1000] = {0};
    	for (int j=0; j<1000; j++){
        	spectra_energy[j] = intensity_input[j][0];
        	spectra_intensity[j] = scaling_factor * intensity_input[j][1]; }

    	gsl_interp_accel *accel_ptr_inten;
    	gsl_spline *spline_ptr_inten;
    	accel_ptr_inten = gsl_interp_accel_alloc ();
    	spline_ptr_inten = gsl_spline_alloc (gsl_interp_cspline, 1000);
    	gsl_spline_init (spline_ptr_inten, spectra_energy, spectra_intensity, 1000);

	//This loop switches the dust composition used in the calculations.
	//It is presently set to calculate all values by default
	for (int x = 0; x < 3; x++ ){

		int i_hi = 0;
 		int i_lo = 0;
  		int proc_number;
  		int proc_remain;
  		int task_proc;
  		int task_remain;

  		int task_number = 91; 
  		MPI_Comm_size( MPI_COMM_WORLD ,&proc_number );
            
  		int task_array[proc_number][3];

  		task_remain = task_number;
  		proc_remain = proc_number;

  		for ( int proc = 1; proc <= proc_number; proc++ ){
    			task_proc = i4_div_rounded ( task_remain, proc_remain );

    			proc_remain = proc_remain - 1;
    			task_remain = task_remain - task_proc;

    			i_lo = i_hi + 1;
    			i_hi = i_hi + task_proc;
    
			task_array[proc-1][0] = proc-1;
			task_array[proc-1][1] = i_lo - 1;
			task_array[proc-1][2] = i_hi - 1;
    		}
		//////////////////////////////////

    		for (int b=0; b < proc_number; b++){
			int order = task_array[b][0]; 
	
			if ( rank == order ){
		 		i_lo = task_array[b][1]; 
		 		i_hi = task_array[b][2];
			}
    		}

    		while (i_lo <= i_hi ){
    
    			double energy = energy_start;
    				
			double grain_size;
			grain_size = ( i_lo )/10.0 + 1.0;
   			string grain_string = static_cast<ostringstream*>( &(ostringstream() << grain_size) )->str();
				
			string input = "../Inputs/" + element[x] + "_Mie_Data/" + element[x] + 
			"_Mie_" + grain_string + "nm.dat";

			fstream input_file(input.c_str(), fstream::in);
    				double cross_input[272][182] = {0};
    				for( int i=0; i<272; i++ ){
        				for( int j=0; j<182; j++ ){
            					input_file >> cross_input[i][j];
        				}
    				}
    				input_file.close();

    			//Calculates the angle between between the Sun and detector from comet center
    			double theta = acos( (pow(rd,2.0) + pow(rg,2.0) - pow(Ro,2.0) ) / (2*rg*rd) );

			//Defines an output file
    			string final_name = "../Results/" + element[x] + "_Results/" + comet_name + "/" 
			+ element[x] + "_output_" + grain_string + "nm_" + comet_name + ".dat";

			int z = 1; 
    			//This loop calculates the intensity emitted for each grain size over the desired energy range
    			while ( energy <= energy_end ){

				double cross_angle[181] = {0};
    				double cross_sec[181] = {0};

    				for ( int k=1; k < 182; k++){
	    				cross_angle[k-1] = cross_input[0][k];

	    				cross_sec[k-1] = cross_input[z][k] / 1e4;	//converts to m^2
    				}			

				//Establishes a interpolation function for the cross-sections
    				gsl_interp_accel *accel_ptr;
    				gsl_spline *spline_ptr;
    				accel_ptr = gsl_interp_accel_alloc ();
    				spline_ptr = gsl_spline_alloc (gsl_interp_cspline, 181);
    				gsl_spline_init (spline_ptr, cross_angle, cross_sec, 181);

				//Scattered differential cross section
        			double scat_cross = gsl_spline_eval (spline_ptr, (180.0/M_PI)*(M_PI-theta), accel_ptr);

        			//Dust particle radius distribution
				double n = abs(alpha - 1) * Np/(a_low) * pow( (a_low/(grain_size*1e-9)), alpha );

        			//Calculates intensity ratio
       				double ratio = n * scat_cross * pow(Ro,2.0)/( pow(rg,2.0)*pow(rd,2.0) ) * delta_a;

				double solar_intensity_temp = gsl_spline_eval (spline_ptr_inten, energy, accel_ptr_inten);
					
				double Intensity = ratio * solar_intensity_temp;

    				//Calculates Brightness of photons scattered per second
    				double Brightness = 0;
    				double diff_cross_Earth = gsl_spline_eval (spline_ptr, (180.0/M_PI)*(M_PI-theta), accel_ptr);

    				for( int l = 0; l<=179; l++ ){
	    				double diff_cross_i = gsl_spline_eval (spline_ptr, l, accel_ptr);

	    				Brightness = Brightness + ( 2*M_PI*pow(rd,2.0)*Intensity )/( diff_cross_Earth ) 
                            * cos((l+0.5)*M_PI/180) * ( sin((l+1)*M_PI/180) - sin(l*M_PI/180) ) * diff_cross_i ;
    				}
                    
                    Brightness = Brightness*(1e4);
                    //Converts Brightness to cm^2 to properly cancel with intensity
                    
    				gsl_spline_free (spline_ptr);
    				gsl_interp_accel_free (accel_ptr);

    				fstream output_file(final_name.c_str(), fstream::out | fstream::app);

    				output_file << scientific << energy << " " << Intensity << " " << Brightness << endl;
    				output_file.close();


    				energy = energy + energy_step;
				z++; 
    			}
			i_lo++;
    		}	

		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0){
   			output_compile(energy_start, energy_end, energy_step, comet_name, x); }

	}

	MPI_Barrier(MPI_COMM_WORLD);
		
	if (rank == 0){
		dust_mixing_ratio(energy_start, energy_end, energy_step, comet_name); }

	gsl_spline_free (spline_ptr_inten);
   	gsl_interp_accel_free (accel_ptr_inten); 

	MPI_Barrier(MPI_COMM_WORLD);
    	return 0;  
}

