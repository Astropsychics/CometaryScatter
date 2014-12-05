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
#include "../dust_large_grain.h"
#include "../dust_input.h"

using namespace std;

int dust_core(double energy_start, double energy_end, double energy_step, int comet_number){

    //inputs number of dust types from comet_input.h
    int dust_types = input_dust_types;

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD , &rank );

    //inputs data from comet_input.h
	double Ro = input_Ro;
	double Q = input_Q[comet_number - 1];
	double vel = input_vel[comet_number - 1];
	double r = input_r[comet_number - 1];
	double rd = input_rd[comet_number - 1];
	double rg = input_rg[comet_number - 1];
	double scaling_factor = input_scaling_factor[comet_number - 1];
	string comet_name = input_comet_name[comet_number - 1];

    //inputs dust model variables from dust_input.h
    double alpha = input_alpha;
    double beta = input_beta;
    double a_low = input_a_low;
    double a_max = input_a_max;
    double delta_a = input_delta_a;
    double gas_mass = input_gas_mass;
    double dust_density = input_dust_density;

	/////////////////////////////////
    //double Np = (Q/4.76e27)*1e27;      //Number of particles, as determined from Fink & Rubin 2012,
	/////////////////////////////////

    //Integrates to find the average particle radius given the selected distribution
    //Note: the equation is a simplified notation of int[a*n(r,a)] / int[n(r,a)]
    double a_avg = 0;
    for( double ap = a_low; ap <= a_max; ap+=a_low ){
        a_avg += (alpha - 1) * pow(a_low/ap, alpha - 1) * a_low; }

    //Calculates dust particle loss rate assuming mass loss rate is equal to
    //the gas mass loss rate
    double Qd = (Q * gas_mass) / ( (4/3) * M_PI * a_avg*a_avg*a_avg * dust_density);

    //Calculates total number of dust particles
    double Np = Qd * r/vel;


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
	for (int x = 0; x < dust_types; x++ ){

        //input parameters for elements from comet_input.h
        string element = input_dust_element[x];

        //determines optimal delegation of jobs to number of cores
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


            //inputs cross section files
			string input = "../Inputs/" + element + "_Mie_Data/" + element +
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
    			double theta = acos( (rd*rd + rg*rg - Ro*Ro ) / (2*rg*rd) );

			//Defines an output file
    			string final_name = "../Results/" + element + "_Results/" + comet_name + "/"
			+ element + "_output_" + grain_string + "nm_" + comet_name + ".dat";

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
       				double ratio = n * scat_cross * Ro*Ro/( rg*rg * rd*rd ) * delta_a;

                    double solar_intensity_temp = gsl_spline_eval (spline_ptr_inten, energy, accel_ptr_inten);

                    double Intensity = ratio * solar_intensity_temp;

    				//Calculates Brightness of photons scattered per second
    				double Brightness = 0;
    				double diff_cross_Earth = gsl_spline_eval (spline_ptr, (180.0/M_PI)*(M_PI-theta), accel_ptr);

    				for( int l = 0; l<=179; l++ ){
	    				double diff_cross_i = gsl_spline_eval (spline_ptr, l, accel_ptr);

	    				Brightness = Brightness + ( 2*M_PI*rd*rd*Intensity )/( diff_cross_Earth )
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

        //Runs output_compile.h
		if (rank == 0){
   			output_compile(energy_start, energy_end, energy_step, comet_name, x); }

	}

	MPI_Barrier(MPI_COMM_WORLD);

    //Pushes results to dust_mixing_ratio.h and dust_large_grain.h
	if (rank == 0){
		dust_mixing_ratio(energy_start, energy_end, energy_step, comet_name);
        dust_large_grain(energy_start, energy_end, energy_step, comet_name);
    }

	gsl_spline_free (spline_ptr_inten);
   	gsl_interp_accel_free (accel_ptr_inten);

	MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}
