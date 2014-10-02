#include <iostream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>

#include "../comet_input.h"
#include "../constants_input.h"
#include "../atomic_cross_section_input.h"
#include "../gas_mixing_ratio.h"

using namespace std;

int gas_core(double energy_start, double energy_end, double energy_step, int comet_number){

    //Reads input from constant_input.h
	double Ro = input_Ro; 

    //Reads inputs from comet_input.h
	double Q = input_Q[comet_number - 1]; 
	double vel = input_vel[comet_number - 1];
	double r = input_r[comet_number - 1];
	double rd = input_rd[comet_number - 1];
	double rg = input_rg[comet_number - 1];
	double scaling_factor = input_scaling_factor[comet_number - 1];
	string comet_name = input_comet_name[comet_number - 1];

    //Calculates the between between the Sun and detector from comet center
    double theta = acos( (pow(rd,2.0) + pow(rg,2.0) - pow(Ro,2.0)) / (2*rg*rd) );

    //Calculates intensity ratio integral
   	 double ratio = Q * (r/vel) * 1/(6*M_PI) * (1 + pow(cos(M_PI - theta),2.0)) * pow(Ro,2.0)/( pow(rg,2.0)*pow(rd,2.0) );

    //Calculates Brightness of photons scattered per second
    //Note: This ratio needs to be multiplied by the intensity ratio to be complete,
    //which it is for each atom in the sections below
    double Brightness = 0.0;
    for( int i = 0; i<=179; i++ ){
    Brightness = Brightness + ( 2*M_PI*pow(rd,2.0) )/( 1 + pow(cos(M_PI - theta), 2.0) )
        * cos((i+0.5)*M_PI/180) * ( sin((i+1)*M_PI/180) - sin(i*M_PI/180) )
		* ( 1 + pow(cos((i+0.5)*M_PI/180),2.0) ) ; }

    Brightness = Brightness*(1e4);
    //Converts Brightness to cm^2 to properly cancel with intensity in next module


    //Pushes data to atomic_cross_section_input.h
    atomic_cross(energy_start, energy_end, energy_step, comet_name, scaling_factor, ratio, Brightness);

    //Pushes data to gas_mixing_ratio.h
	gas_mixing_ratio(energy_start, energy_end, energy_step, comet_name);
		
   	return 0; 
}

//Version history
//
