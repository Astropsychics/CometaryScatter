#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>

#include "atomic.h"
#include "chandra.h"
#include "dust_scatter.h"
#include "fluo.h"
#include "cx_calculations.h"
#include "cx_compare.h"

using namespace std;

//Energy range (in keV)
double energy_start = 0.300;
double energy_end = 3.000;
double energy_step = 0.010;

//energy peak width (in keV)
double width = 0.050;

int main(int argc, char *argv[]){

    MPI::Init();

	int comet_number = 0;
	bool comet_all = false;

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD , &rank );
	if (rank == 0){

		//Prompt user to select which comet to analyze
        cout << "\n" << "Welcome to the Cometary Emission Model. \n \n";
		cout << "Please select the number of which comet you'd like to analyze. \n";
		cout<< "1.8P  2.Encke  3.IZ  4.LS4 \n" << "5.MH  6.ISON  7.PanSTARRS  8.all comets \n";
		cin >> comet_number;
        cout << endl;

		if ( 1 > comet_number || comet_number > 8 ){
			cout << "I'm sorry, Dave. I can't do that. \n";
			return 0; }

        //performs cx calculations that will be scaled to each comet later
        cx_calculations(energy_start, energy_end, energy_step, width);
	}

	MPI_Bcast(&comet_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	if (comet_number == 8 ){
        comet_all = true;
		comet_number = 1; }


	while ( comet_number < 8 ){

        MPI_Barrier(MPI_COMM_WORLD);

        //Pushes data to atomic.h, fluo.h, cx_calculations.h 
        //This is performed by one core as calculations are trivial
        if (rank == 0){
        		gas_core(energy_start, energy_end, energy_step, comet_number);
        		fluo_core(energy_start, energy_step, comet_number);
                cx_compare(energy_start, energy_end, energy_step, comet_number); }

        //Pushes data to dust_scatter.h
		dust_core(energy_start, energy_end, energy_step, comet_number);

        //Pushes data to chandra.h
        //Calculations performed by one core, again, due to triviality
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0){
            chandra_core(energy_start, energy_end, energy_step, width, comet_number); }

        //Add +1 to counter to begin next comet, else ends program
        if (comet_all == true){
            if (rank ==0) cout << comet_number << " of 7 completed. \n";
            comet_number++; }
        else {
			MPI_Barrier(MPI_COMM_WORLD);
			MPI::Finalize();
            return 0; }
	}

    MPI_Barrier(MPI_COMM_WORLD);
	MPI::Finalize();
   	return 0;
}
