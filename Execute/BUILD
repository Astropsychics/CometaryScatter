#!/bin/bash

mkdir ../Results/C_Results
mkdir ../Results/H2O_Results
mkdir ../Results/Si_Results

mkdir ../Results/C_Results/8P
mkdir ../Results/C_Results/Encke
mkdir ../Results/C_Results/IZ
mkdir ../Results/C_Results/LS4
mkdir ../Results/C_Results/MH
mkdir ../Results/C_Results/ISON
mkdir ../Results/C_Results/PanSTARRS

mkdir ../Results/H2O_Results/8P
mkdir ../Results/H2O_Results/Encke
mkdir ../Results/H2O_Results/IZ
mkdir ../Results/H2O_Results/LS4
mkdir ../Results/H2O_Results/MH
mkdir ../Results/H2O_Results/ISON
mkdir ../Results/H2O_Results/PanSTARRS

mkdir ../Results/Si_Results/8P
mkdir ../Results/Si_Results/Encke
mkdir ../Results/Si_Results/IZ
mkdir ../Results/Si_Results/LS4
mkdir ../Results/Si_Results/MH
mkdir ../Results/Si_Results/ISON
mkdir ../Results/Si_Results/PanSTARRS

make clobber
make
make clean 
mpirun -np $1 ./MAIN

rm -rf ../Results/C_Results
rm -rf ../Results/H2O_Results
rm -rf ../Results/Si_Results
