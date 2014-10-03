#ifndef COMET_INPUT_H
#define COMET_INPUT_H

using namespace std;

// array[0] = 8P
// array[1] = Encke
// array[2] = IZ
// array[3] = LS4 
// array[4] = MH
// array{5] = ISON
// array[6] = PanStarrs

//comet name
string input_comet_name[7] = {"8P","Encke","IZ","LS4","MH","ISON","PanSTARRS"};

//production rate, in mol/sec
double input_Q[7] = {2.2e28, 0.7e28, (1.9e29)/pow(0.81,4.0), 3e28, 2.0e29, 1e29, 1e29};

//outgas velocity, in m/sec
double input_vel[7] = {800, 800, 800, 800, 800, 800, 800};

//atmosphere radius, in m
double input_r[7] = {1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9};

//comet-Earth distance, in m
double input_rd[7] = {0.25*1.496e11, 0.28*1.496e11, 0.45*1.496e11, 0.53*1.496e11, 1.37*1.496e11,
    0.95*1.496e11, 1.442*1.496e11};

//Sun-comet distance, in m
double input_rg[7] = {1.10*1.496e11, 0.89*1.496e11, 0.81*1.496e11, 0.80*1.496e11, 1.26*1.496e11,
   1.176*1.496e11, 1.102*1.496e11};

//solar spectrum scaling factor, based on GOES observations
double input_scaling_factor[7] = {(0.134 * 6.24e4)/2.36e5, (0.134 * 4.37e5)/2.36e5,
        (0.134 * 3.00e6)/2.36e5, (0.134 * 3.12e6)/2.36e5, (0.134 * 6.24e5)/2.36e5,
        (0.134 * 1.20e6)/2.36e5, (0.134 * 3.75e5)/2.36e5};

//exposure time, in sec
double input_t_exp[7] = {47000, 44000, 24000, 9400, 16900, 27000, 30360};

////////////////////////////////////////

//number of different gases included
int input_gas_types = 4;

//input string of different gases included
string input_gas_element[4] = {"H","C","N","O"};

//atomic mixing ratio
double input_gas_mixing_ratio[4] = {2, 0.1, 0.05, 1};

///////////////////////////////////////////

//number of different dusts included
int input_dust_types = 3;

//input string of different dust types included
string input_dust_element[3] = {"C","H2O","Si"};
string input_dust_element_revised[3] = {"C_dust","H2O","Si_dust"};

//dust mixing ratio
double input_dust_mixing_ratio[3] = {0.1, 0.85, 0.05};

#endif