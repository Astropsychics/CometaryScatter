#include <iostream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>

#include "../comet_input.h"
#include "../constants_input.h"
#include "../fluorescence_calculations.h"

using namespace std;

int fluo_core(double energy_start, double energy_step, int comet_number){

	double Ro = input_Ro;
	double Q = input_Q[comet_number - 1]; 
	double vel = input_vel[comet_number - 1];
	double r = input_r[comet_number - 1];
	double rd = input_rd[comet_number - 1];
	double rg = input_rg[comet_number - 1];
        double scaling_factor = input_scaling_factor[comet_number- 1];
	string comet_name = input_comet_name[comet_number - 1];


        //Calculates intensity ratio integral
        double fluo_ratio = Q * (r/vel) * 1/(4*M_PI) * pow(Ro,2.0)/( pow(rg,2.0)*pow(rd,2.0) );

        fluorescence_calc(energy_start, energy_step, comet_name, scaling_factor, fluo_ratio);

	return 0; 
}


//Version history
// 
