 #include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>

#include "../cx_input.h"

using namespace std;

int cx_calculations(double energy_start, double energy_end, double energy_step, double width){

    //calculates numbers of steps and creates appropriately-sized array
    double energy = energy_start;
    int energy_row = (energy_end - energy_start)/energy_step + 1;


    vector<double> input_energy(energy_row, 0);

    for ( int x = 0; x < energy_row; x++ ){
        input_energy[x] = energy;
        energy += energy_step; }

    //input effective area
    fstream input_area_file("../Inputs/Effective_Area/acis_effective_area.dat", fstream::in);

    int area_row;
    input_area_file >> area_row;

    double input_area[area_row][3];
    for( int i=0; i<area_row; i++ ){
        for( int j=0; j<3; j++ ){
            input_area_file >> input_area[i][j]; } }
	input_area_file.close();

    //inputs number of elements from cx_input.h
    int cx_elements = input_cx_elements;
    //calculates sum of ratio required for future normalization
    double sum = 0;
    for (int a = 0; a < cx_elements; a++) sum += input_cx_sigma[a]*input_cx_n[a];

    //defines cx spectrum array
    vector<double> cx_spectrum(energy_row, 0);

    //calculates cx emissions for each element
    for ( int x = 0; x < cx_elements; x++ ){

        //input strings dependent on whether calculations are performed for high or low solar wind activity
        string sw_hilow;
        string line_name;
        double eta = 0;

        //defines solar wind activity: low = 1 , high = 2
        int sw_activity = 1;

        if ( sw_activity == 1 ){
            line_name = input_cx_line_name_low[x];
            eta = 1/1780.0 * input_cx_n[x] * input_cx_sigma[x] / sum;
            sw_hilow = "Low"; }

        if (sw_activity == 2 ){
            line_name = input_cx_line_name_high[x];
            eta = 1/1550.0 * input_cx_n[x] * input_cx_sigma[x] / sum;
            sw_hilow = "High"; }

        //reads in line spectrum
        string line_string = "../Inputs/CX_Tables/" + sw_hilow + "/" + line_name + ".dat";
        fstream input(line_string.c_str(), fstream::in);

        int row = 0;
        input >> row;

        double input_cx[row][2];
        for ( int i=0; i<row; i++ ){
            for( int j=0; j<2; j++ ){
                input >> input_cx[i][j]; }

            input_cx[i][0] = input_cx[i][0] / 1000.0;
            input_cx[i][1] = eta * input_cx[i][1];}
        input.close();

        //creates gaussian functions for each cx line
        for (int p = 0; p < row; p++){
            for ( int q = 0; q < energy_row; q++ ){

                double peak = input_cx[p][0];
                double area = input_cx[p][1];
                double amplitude = area / ( width * sqrt(2*M_PI) );

                //define a Gaussian function for the cx peaks
                cx_spectrum[q] += amplitude * exp( -pow(input_energy[q] - peak,2.0) / (2 * width*width) );
            }
        }
    }

    //declares output file names
    ofstream spectra("../Results/CX_spectrum.dat");
    ofstream Chandra("../Results/CX_Chandra_spectrum.dat");

    for (int l=0; l<energy_row; l++){
        int z = 0;
        while ( input_area[z][1] <= input_energy[l] ) z++;

        //outputs total intensity spectrum to a file
        spectra << scientific << input_energy[l] << " " << cx_spectrum[l] << endl;

        Chandra << scientific << input_energy[l] << " " << cx_spectrum[l] * input_area[z][2] << endl;

    }
    spectra.close();
    Chandra.close();

    return 0;
}
