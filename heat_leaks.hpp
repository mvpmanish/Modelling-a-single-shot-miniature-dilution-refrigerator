//Header file to model the various heat leaks from the fridge i.e. the number of Joules per second that are lost in the operation of the fridge.

#ifndef HEAT_LEAKS_HPP
#define HEAT_LEAKS_HPP

#include <cmath>

#include "mixing_chamber.hpp"

using namespace std;

//---------------------------------------------------RADIATIVE HEAT LEAKS-------------------------------------------------
//Constants
extern double e_cu, e_st;  //Emissivity of copper (assume 0.7) and steel (assume 0.5)
extern double sb_const;  //Stefan-Boltzmann constant [Wm^-2K^-4]
extern double still_temp;  //Teperature of the still

//Radiative heat transfer from the still to the mixing chamber
double still_to_mc_heat_leak(double &T);

//Radiative heat transfer from the bottom of the cryostat to the mixing chamber
double cryo_to_mc_heat_leak(double &T);

//--------------------------------------------------CONDUCTIVE HEAT LEAKS------------------------------------------------
//Constants
extern double cf_leg_area;  //Area of one carbon-fibre support leg [m^2]
extern double p1_cf, p2_cf;  //Constants which are from linear fit of heat leak vs temperature difference p1 in [Wm^-2K^-1] and p2 in [Wm^-2]
extern double pi;

//Conductive heat leak from carbon fiber hexapod legs, T in [K] and returns in [W]
double cf_heat_leak(double &T);  

//Heat leak through the steel capillary, T in [K] and returns in [W]
double steel_capillary_heat_leak(double &T);  

//Heat leak from the superfluid Helium-II in the capillary
double helium_capillary_heat_leak(double &T);  

//---------------------------------------OVERALL HEAT LEAK---------------------------------------------------------------
//Total heat leak in [W]
double q_leak(double &T);

//---------------------------------------------HE-II PROPERTIES NEEDED----------------------------------------------------
//Viscosity of He-II [Pa*s]
double He_viscosity(double &T);

//Entropy in the dilute phase 
double dil_He_entropy(double &T);


#endif HEAT_LEAKS_HPP