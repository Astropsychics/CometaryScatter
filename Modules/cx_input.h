#ifndef CX_INPUT_H
#define CX_INPUT_H

using namespace std;

//number of elements
int input_cx_elements = 20;

//element order: C6+, C5+, N7+, N6+, N5+, O8+, O7+, O6+, Ne9+, Ne8+, Mg10+, Mg9+,
//               Fe13+, Fe12+, Fe11+, Fe10+, Si9+, S10+, S9+, S8+

//cross section, in units of 10-15 cm^2
double input_cx_sigma[20] = {4.16, 0.722, 9.5, 3.71, 2.27, 6.3, 5.3, 3.67, 7.2, 3.7, 3.7,
                            4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0};

//particle density, relative to O
//double input_cx_n[20] = {0.318, 0.210, 0.006, 0.058, 0.065, 0.070, 0.200, 0.730, 0.004,
//                         0.084, 0.098, 0.052, 0.002, 0.007, 0.023, 0.031, 0.021, 0.005, 0.016, 0.019};


//Increased O8+ production
//double input_cx_n[20] = {0.318, 0.210, 0.006, 0.058, 0.065, 0.130, 0.200, 0.670, 0.004,
//                         0.084, 0.098, 0.052, 0.002, 0.007, 0.023, 0.031, 0.021, 0.005, 0.016, 0.019};

//Increased O8+ and Ne9+ production
double input_cx_n[20] = {0.318, 0.210, 0.006, 0.058, 0.065, 0.130, 0.200, 0.670, 0.020,
                         0.084, 0.098, 0.052, 0.002, 0.007, 0.023, 0.031, 0.021, 0.005, 0.016, 0.019};


//cx low sw speed table filename
string input_cx_line_name_low[20] = {"C6collisionwithH2V0.19","C5collisionwithH2V0.18",
    "N7collisionwithH2OV0.20","N6collisionwithH2V0.21","N5collisionwithH2V0.19",
    "myGreenO8H2O","O7collisionwithH2OgreenfitV0.20","O6collisionwithHV0.19",
    "Ne9collisionwithH2V0.20","Ne8CollisionwithH2V0.25",
    "Mg10CollisionwithH2V0.25","Mg9collisionwithH2","Fe13collisionwithH2",
    "Fe12collisionwithH2","Fe11collisionwithH2","Fe10collisionwithH2",
    "Si9collisionwithH2","S10collisionwithH2","S9collisionwithH2",
    "S8collisionwithH2"};

//cx high sw speed table filename
string input_cx_line_name_high[20] = {"xxxxxxxxx","xxxxxxxxx","xxxxxxxxx",
    "xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx",
    "xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx",
    "xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx"};

#endif
