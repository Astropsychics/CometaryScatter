#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <vector>

#include "../comet_input.h"

using namespace std;

//Energy width (in keV)
double width = 0.050;

int chandra_core(double energy_start, double energy_end, double energy_step, int comet_number){
    
    	//input effective area
    	fstream input_area_file("../Inputs/Effective_Area/effective_area_Chandra.dat", fstream::in);
    
    	int area_row;
    	input_area_file >> area_row;
    
    	double input_area[area_row][3];
    	for( int i=0; i<area_row; i++ ){
        	for( int j=0; j<3; j++ ){
            		input_area_file >> input_area[i][j]; } }
	input_area_file.close();
        
	double t_exp = input_t_exp[comet_number - 1];
	string comet_name = input_comet_name[comet_number - 1];
        
	int row = (energy_end - energy_start)/energy_step + 1;
	vector<double> total_inten(row, 0);
        
	double input_inten[row][3];
        double input_energy[row];
        
        string element_type[2] = {"gas","dust"};
        for( int x = 0; x < 2; x++ ){
            
		string Inten = "../Results/" + comet_name + "/" + element_type[x] +
                "_total_output_" + comet_name + ".dat";
		fstream input(Inten.c_str(), fstream::in);
            
            	for( int i=0; i<row; i++ ){
                	for( int j=0; j<3; j++ ){
                    	input >> input_inten[i][j];} }
            	input.close();
        
            	for ( int y = 0; y < row; y++ ){
                	input_energy[y] = input_inten[y][0];
                	total_inten[y] = total_inten[y] + input_inten[y][1]; }
        }

        string Inten_fluo = "../Results/" + comet_name + "/fluorescence_total_output_" + comet_name + ".dat";
        fstream input_fluo(Inten_fluo.c_str(), fstream::in);
        
        int fluo_row;
        input_fluo >> fluo_row;
        
        double input_intenfluo[fluo_row][3];
        for ( int i=0; i<fluo_row; i++ ){
            	for( int j=0; j<3; j++ ){
            		input_fluo >> input_intenfluo[i][j]; } }
        input_fluo.close();
        
	vector<double> fluo_spectrum(row, 0);
	for (int p=0; p<fluo_row; p++){
           	for (int q=0; q<row; q++){

	   		double peak = input_intenfluo[p][1];
           		double area = input_intenfluo[p][2];
           		double amplitude = area / ( width * sqrt(2*M_PI) );
  	
           		//define a Gaussian function for the fluorescence peaks then add it with the scattering components
           		fluo_spectrum[q] = fluo_spectrum[q] + amplitude 
			* exp( -pow(input_energy[q] - peak,2.0) / (2 * pow(width,2.0)) );
           	}
        }

        string spectra_name = "../Results/" + comet_name + "/spectra_total_" + comet_name + ".dat";
        ofstream spectra(spectra_name.c_str());
        
        string Chandra_name = "../Results/" + comet_name + "/Chandra_spectra_total_" + comet_name + ".dat";
        ofstream Chandra(Chandra_name.c_str());
        
        for (int l=0; l<row; l++){
            	int z = 0;
            	while ( input_area[z][1] <= input_energy[l] ){
                	z++; }

            	total_inten[l] = total_inten[l] + fluo_spectrum[l];

            	//outputs total intensity spectrum to a file
            	spectra << scientific << input_energy[l] << " " << total_inten[l] << endl;
            
           	 Chandra << scientific << input_energy[l] << " " << total_inten[l] * input_area[z][2] * t_exp << endl;
        }
        spectra.close();
        Chandra.close();

    	return 0; 
}
