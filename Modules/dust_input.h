#ifndef DUST_INPUT_H
#define DUST_INPUT_H

using namespace std;

double input_alpha = 2.5;           //Dust density dependence on grain radius
double input_beta = 2.0;            //Dust density dependence on atmosphere radius
double input_gamma = 2.0;           //Dust density dependence on cross section

double input_a_low = 1e-9;		    //Lower bound of grain size integration
double input_delta_a = 0.1e-9;	    //Grain radius step size

double input_a_max = 1e-2;          //Maximum grain size, used for large grain calculations

#endif
