#include "mpi.h"
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <stdlib.h>

#include "atomic.h" 
#include "chandra.h" 
#include "dust_scatter.h"
#include "fluo.h" 

using namespace std;

//Energy range (in keV)
double energy_start = 0.300;
double energy_end = 3.000; 
double energy_step = 0.010; 

int main(int argc, char *argv[]){
    
    MPI::Init();
    
	int comet_number = 0;
	bool comet_all = false;
	
	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD , &rank );
	if (rank == 0){
        
		//Prompt user to select which comet to analyze
		cout << "Please select the number of which comet you'd like to analyze." << endl;
		cout<< "1.8P  2.Encke  3.IZ  4.LS4  5.MH  6.all comets" << endl;
		cin >> comet_number;
        
		if ( 1 > comet_number || comet_number > 6 ){
			cout << "I'm sorry, Dave. I can't do that." << endl;
			return 0; }
	}
	
	MPI_Bcast(&comet_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
    
	if (comet_number == 6 ){
        comet_all = true;
		comet_number = 1; }
    
    
	while ( comet_number < 6 ){

        	MPI_Barrier(MPI_COMM_WORLD);
		if (rank == 0){
        		gas_core(energy_start, energy_end, energy_step, comet_number);
        		fluo_core(energy_start, energy_step, comet_number); }
       
		dust_core(energy_start, energy_end, energy_step, comet_number);
        
        	MPI_Barrier(MPI_COMM_WORLD);
        	if (rank == 0){
        		chandra_core(energy_start, energy_end, energy_step, comet_number); }  

        	if (comet_all == true){
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
