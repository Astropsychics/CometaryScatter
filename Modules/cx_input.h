#ifndef CX_INPUT_H
#define CX_INPUT_H

using namespace std;

//number of elements
int input_cx_elements = 11;

//cross section, in units of 10-17 cm^2
double input_cx_sigma[11] = {4.16, 0.722, 9.5, 3.71, 2.27, 6.3, 5.3, 3.67, 7.2, 3.7};

//particle density, relative to O
double input_cx_n[11] = {.08, .20, .017, .035, .065, .014, .08, .73, .002, .094, .098};

//cx low sw speed table filename
string input_cx_line_name_low[11] = {"C6collisionwithH2V0.19","C5collisionwithH2V0.18",
    "N7collisionwithH2OV0.20","N6collisionwithH2V0.21","N5collisionwithH2V0.19",
    "myGreenO8H2O","O7collisionwithH2OgreenfitV0.20","O6collisionwithHV0.19",
    "Ne9collisionwithH2V0.20","Ne8CollisionwithH2V0.25",
    "Mg10CollisionwithH2V0.25"};
    
//cx high sw speed table filename
string input_cx_line_name_high[11] = {"xxxxxxxxx","xxxxxxxxx","xxxxxxxxx",
    "xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx","xxxxxxxxx",
    "xxxxxxxxx","xxxxxxxxx","xxxxxxxxx"};

#endif
