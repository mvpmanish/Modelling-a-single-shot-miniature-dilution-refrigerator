//header file to define the vapour pressure curves of helium

#ifndef VAPOUR_PRESSURE_HPP
#define VAPOUR_PRESSURE_HPP

#include<vector>
#include<cmath>
#include<string>
#include<fstream>
#include<iostream>

using namespace std;

//------------------------------------------------FUNCTIONS---------------------------------------------
extern vector <double> vp_He3_coeffs; // units of [mBar]

//Function to calculate vapour pressure in [Pa]
double vap_pressure_He3(double T);

//Function to calculate vapour pressure of Helium-4 in [Pa] from [Boreghi 1998]
double vap_pressure_He4(double T);

//Function to write vapour pressure vs Temperature to a file
void write_vap_pressure_to_file(string filename, double start_temp, double end_temp, double temp_interval);

#endif