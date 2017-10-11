// Header file to define all the things happening in the still
// Manny & Chief

#ifndef STILL_HPP
#define STILL_HPP

#define PI 3.14159265359

#include<vector>
#include<cmath>

#include "vapour_pressure.hpp"

using namespace std;

//-----------------------------------------------------------CONSTANTS-----------------------------------------------------------------------
double molar_mass_He3{ 3.016e-3 };  //[kg/mol]
double length_of_pipe{ 175e-3 };  //[m]

//----------------------------------------------------------CLASS--------------------------------------------------------------
class still
{
private:
	double T; // Temperature  [K]
	double C; // Heat Capacity  [J/K]
	double V; // Volume [m^3]
	double P; // Pressure at top of still [Pa]
	double diameter{ 6.05e-3 };  //units of [m]
	double n3_dot; //He-3 molar flow rate [mol/s]
	
public:
	still() {}; // Default constructor
	still(double Ti, double V); // Param constructor
	~still() {}; // Default destructor

	double molecular_flow(double T, double pressure, double molar_mass_in, double diameter_in, double length_in);  //Molecular flow rate [mol/s]

	// Access functions
	double temp() { return T; }
	double specfic_heat() { return C; }
	double vol() { return V; }
	double get_n3_dot() { return n3_dot; }
	
	void set_n3_dot();  //Function to set n3_dot for the still
};

//Calculates the moelcular flow rate [mol/s] given certain parameters
double still::molecular_flow(double Ti, double vap_pressure, double molar_mass, double diameter, double length)
{
	double r = diameter/2;
	double n_dot
	{
		(4 / 3)*sqrt((2 * PI) / (8.31*Ti*molar_mass)) * (pow(r,3) / length) * vap_pressure_He3(Ti)
	};
	return n_dot;
}

//Parameterised constructor
still::still(double Ti, double Vi)
{
	if (Ti < 0.28)  //Vapour pressure for He3 not valid below this temperature
	{
		cerr << endl << "Still temperature is too low for model! ( < 0.28K)" << endl;
		return;
	}
	T = Ti; V = Vi;
	P = vap_pressure_He3(Ti);
	n3_dot = molecular_flow(Ti, P, molar_mass_He3, diameter, length_of_pipe);
}

//Function to set the value of He-3 molar flow rate [mol/s] for the still
void still::set_n3_dot()
{
	n3_dot = molecular_flow(T, P, molar_mass_He3, diameter, 175e-3);  //175e-3 = length of the tube [m]
}


#endif