//This is a header file to define a class for the mixing chamber in the single-shot miniature dilution fridge model
//Manish Patel & Adriano Parisi
//Last date modified: 25/10/2016

#ifndef MIXING_CHAMBER_HPP
#define MIXING_CHAMBER_HPP

#include<vector>

#include "vapour_pressure.hpp"

using namespace std;

//---------------------------------------------------------CONSTANTS-------------------------------------------------------

//Declaration of external constants 
extern double x_t, T_t;  //Tricritical concentration and temperature
extern double h_at_x_t, dh_dx;
extern double He3_mass, He4_mass;  //Atomic masses of He-3 and He-4 [kg]
extern double N_A;  //Avagadro's number [/mol]
extern double gamma_e;  //Electronic heat capacity coefficient

//Vectors of constants for fit determining specific heat of helium-3/helium-4 in two phase region [Chaudhry 2009]
extern vector<double> A2_1;  //T < 0.15K
extern vector<double> A2_2;  // T > 0.15K
extern vector<double> B2_1;  //T < 0.15K
extern vector<double> B2_2;  //T >0.

//--------------------------------------------TWO-PHASE FUNCTIONS---------------------------------------------------
//Dilute phase He-3 concentration
double two_phase_x_dilute(double &T); 

//Concentrated phase He-3 concentration
double two_phase_x_conc(double &T);

//Specific heat in the two phase region [J/(K*kg)]
double two_phase_specific_heat(double &T);

//----------------------------------------------------CLASS--------------------------------------------------------------

//Mixing chamber class
class mixing_chamber
{
private:
	double T;  //Temperature of chamber in K
	double x_dil, x_conc;  //concentration of He-3 in dilute and concentrated phases
	double C;  //specific heat in J/k
	double n3, n4;  //Number of mols of helium-3 and helium-4
	double V;  //Volume of mixing chamber
	double P;  //Pressure inside mixing chamber
public:
	//Default constructor
	mixing_chamber(){};

	//Parameterised constructor
	mixing_chamber(double Ti, double Vi, double n3i, double n4i);

	//Default destructor
	~mixing_chamber(){};

	//Access functions
	double temp(){ return T; }
	double dilute_conc(){ return x_dil; }
	double conc_conc(){ return x_conc; }
	double specific_heat(){ return C; }
	double mol_He3(){ return n3; }
	double mol_He4(){ return n4; }
	double vol(){ return V; }
	double pressure(){ return P; }

	//Enthalpy
	double h_dil();  //Enthalpy of dilute phase
	double h_conc(); //Enthalpy of concentrated phase
	double delta_h(){ return h_conc() - h_dil(); }  //Difference in enthalpy between the two phases

	//Function to lower the energy [J] of the mixing chamber 
	void lower_energy(double &E)  //C_copper + mass_mixture*C_mixture = dQ/dT, C_copper = gamma_e*T
	{
		mixing_chamber mc((T - E / (C*(n3*He3_mass*N_A + n4*He4_mass*N_A)+ gamma_e*T)), V, n3, n4);
		T = mc.T;
		x_dil = mc.x_dil; x_conc = mc.x_conc;
		C = mc.C;
		P = mc.P;
	}

	//Function to decrease the mols of He-3 in the mixing chamber
	void lower_n3(double &n3_out)
	{
		n3 = n3 - n3_out;
	}

	//Function to increase number of mols of He-4 in th mixing chamber
	void increase_n4(double &n4_in)
	{
		n4 = n4 + n4_in;
	}

};


#endif