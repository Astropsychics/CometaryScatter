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

    //inputs effective area
    fstream input_area_file("../Inputs/Effective_Area/effective_area_Chandra.dat", fstream::in);

    int area_row;
    input_area_file >> area_row;

    double input_area[area_row][3];
    for( int i=0; i<area_row; i++ ){
        for( int j=0; j<3; j++ ){
            input_area_file >> input_area[i][j]; } }
	input_area_file.close();

    //Inputs data from comet_input.h
	double t_exp = input_t_exp[comet_number - 1];
	string comet_name = input_comet_name[comet_number - 1];

    //Calculates numbers of steps to perform and generates energy and intensity arrays based
    //on the results
	int row = (energy_end - energy_start)/energy_step + 1;
	vector<double> total_inten(row, 0);
    vector<double> dust_inten(row, 0);
    vector<double> gas_inten(row, 0);

	double input_inten[row][3];
    double input_energy[row];

    //inputs gas and dust data
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
            total_inten[y] = total_inten[y] + input_inten[y][1];

            if ( x == 0 ){
                gas_inten[y] = input_inten[y][1];}
            if ( x == 1 ){
                dust_inten[y] = input_inten[y][1];}
        }
    }

    //inputs fluoresence data
    string Inten_fluo = "../Results/" + comet_name + "/fluorescence_total_output_" + comet_name + ".dat";
    fstream input_fluo(Inten_fluo.c_str(), fstream::in);

    int fluo_row;
    input_fluo >> fluo_row;

    double input_intenfluo[fluo_row][3];
    for ( int i=0; i<fluo_row; i++ ){
        for( int j=0; j<3; j++ ){
            input_fluo >> input_intenfluo[i][j]; } }
    input_fluo.close();

    //creates fluorescence spectrum by generating and summing gaussians curves for all peaks
	vector<double> fluo_spectrum(row, 0);
	for (int p=0; p<fluo_row; p++){
        for (int q=0; q<row; q++){

            double peak = input_intenfluo[p][1];
            double area = input_intenfluo[p][2];
            double amplitude = area / ( width * sqrt(2*M_PI) );

            //define a Gaussian function for the fluorescence peaks
            fluo_spectrum[q] = fluo_spectrum[q] + amplitude
            * exp( -(input_energy[q] - peak)*(input_energy[q] - peak) / (2 * width*width) );
        }
    }

    //Creates output files
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

        Chandra << scientific << input_energy[l] << " " << total_inten[l] * input_area[z][2] << endl;
    }
    spectra.close();
    Chandra.close();

    //Outputs gas and dust spectra separately, if required
    int ADDITIONAL_DATA = 1;
    if ( ADDITIONAL_DATA == 1 ){
        string dust_out_name = "../Results/" + comet_name + "/Chandra_dust_spectra_" + comet_name + ".dat";
        ofstream dust_out(dust_out_name.c_str());

        string gas_out_name = "../Results/" + comet_name + "/Chandra_gas_spectra_" + comet_name + ".dat";
        ofstream gas_out(gas_out_name.c_str());

        for (int l=0; l<row; l++){
            int z = 0;
            while ( input_area[z][1] <= input_energy[l] ){
                z++; }

            //outputs intensity spectra to file
            dust_out << scientific << input_energy[l] << " " << dust_inten[l] * input_area[z][2] << endl;

            gas_out << scientific << input_energy[l] << " " << gas_inten[l] * input_area[z][2] << endl;
        }
    dust_out.close();
    gas_out.close();
    }

    return 0;
}
