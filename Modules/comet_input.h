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

string input_comet_name[7] = {"8P","Encke","IZ","LS4","MH","ISON","PanSTARRS"};

double input_Q[7] = {2.2e28, 0.7e28, (1.9e29)/pow(0.81,4.0), 3e28, 2.0e29, 1e29, 1e29};

double input_vel[7] = {800, 800, 800, 800, 800, 800, 800};

double input_r[7] = {1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1e9};

double input_rd[7] = {0.25*1.496e11, 0.28*1.496e11, 0.45*1.496e11, 0.53*1.496e11, 1.37*1.496e11,
    0.95*1.496e11, 1.442*1.496e11};

double input_rg[7] = {1.10*1.496e11, 0.89*1.496e11, 0.81*1.496e11, 0.80*1.496e11, 1.26*1.496e11,
    1.176*1.496e11, 1.102*1.496e11};

double input_scaling_factor[7] = {(0.134 * 6.24e4)/2.36e5, (0.134 * 4.37e5)/2.36e5,
        (0.134 * 3.00e6)/2.36e5, (0.134 * 3.12e6)/2.36e5, (0.134 * 6.24e5)/2.36e5,
        (0.134 * 1.20e6)/2.36e5, (0.134 * 3.75e5)/2.36e5};

double input_t_exp[7] = {47000, 44000, 24000, 9400, 16900, 27000, 30360};

#endif
