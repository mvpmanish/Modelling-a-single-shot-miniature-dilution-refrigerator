//cpp file to define teh functions declared in "heat_leaks.hpp"

#include "heat_leaks.hpp"

//---------------------------------------------------RADIATIVE HEAT LEAKS-------------------------------------------------
//Radiative heat transfer from the bottom of the cryostat to the mixing chamber in [W]
double cryo_to_mc_heat_leak(double &T)
{
	double T_cryo{ 4 };  //Cryostat at 4K
	double q{ (e_cu*e_cu*sb_const*pow(T_cryo, 4) - e_cu*e_cu*sb_const*pow(T, 4)) / (e_cu + e_cu - e_cu*e_cu) };
	double area_bottom_hexapod{0.10521};  //[m^2]
	return q*area_bottom_hexapod;
}

//Radiative heat transfer from the still to the mixing chamber in [W]
double still_to_mc_heat_leak(double &T)
{
	double q{ (sb_const*(pow(still_temp, 4) - pow(T, 4)) / ((1 / e_cu) + (1 / e_cu) - 1)) };
	double area_outer_mc_wall{ 4.04637e-3 };  //in [m^2]
	return q*area_outer_mc_wall;
}

//--------------------------------------------------CONDUCTIVE HEAT LEAKS-------------------------------------------------
double cf_heat_leak(double &T)
{
		double Tend{ 0.7 };  //End temperature of cf leg [K]
		double av_T{ (Tend + T)/2 };  //Average temperatur of leg
		return (p1_cf*av_T + p2_cf)*cf_leg_area * 6;  //Multiply by 6 as there are 6 legs
}

//CHECK THIS AGAIN: GIVES 1mW heat leak! -> too high
double steel_capillary_heat_leak(double &T)
{
	double Tend{ 0.7 };  //End temperature of capillary [K]
	double av_T{ (Tend - T)/2 }; //Average temperature of steel capillary
	double l{ 542e-3 };  //length of capillary
	double k{ 1.69e-3 + 1.12e-1*T };  //Thermal conductivity Wm^-1K^-1
	return k*av_T*l;
}

//Heat leak from the superfluid Helium-II in the capillary
double helium_capillary_heat_leak(double &T)
{
	double T_end{ 0.9 };  //End temperature of capillary [K]
	double l{ 542e-3 };  //length of capillary [m]
	double R{ 0.3e-3 };  //Radius of capillary [m]
	double S{ dil_He_entropy(T) };  //Entropy of helium in capillary [J/K]
	double n{ He_viscosity(T) };  //Viscosity of helium in capillary [Pas]
	double Tm{ (T + T_end) / 2 }, Delta_T{ T_end - T };  //Mean temperature and temperature difference [K]

	double Q{ pi*pow(R, 4)*pow(S, 2)*Tm*Delta_T / (8 * n*l) };
	return Q;

}
//-------------------------------------------OVERALL HEAT LEAK----------------------------------------------------
//Function to sum all individual heat leaks in [W]
double q_leak(double &T)
{
	double q{ cf_heat_leak(T) + /*steel_capillary_heat_leak(T) */+ helium_capillary_heat_leak(T) + cryo_to_mc_heat_leak(T) + still_to_mc_heat_leak(T) + cryo_to_mc_heat_leak(T) };
	return q;
}

//--------------------------------------------HE-II PROPERTIES NEEDED----------------------------------------------
double He_viscosity(double &T)
{
	double n{ exp(-8.8*T - 4) };  //Formula extrapolated from Berenghi
	return n;
}

//Entropy in the dilute phase 
double dil_He_entropy(double &T)
{
	double x_term{ two_phase_x_dilute(T) - x_t };

	double B_term{ 0 };
	for (int i{ 0 }; i < 6; i++)
	{
		if (T < 0.15)
		{
			B_term += B2_1[i] * pow(T, i + 1) / (i + 1);
		}
		else
		{
			B_term += B2_2[i] * pow(T, i + 1) / (i + 1);
		}
	}

	double A_term{ 0 };
	for (int j{ 0 }; j < 6; j++)
	{
		if (T < 0.15)
		{
			A_term += A2_1[j] * pow(T, j + 1) / (j + 1);
		}
		else
		{
			A_term += A2_2[j] * pow(T, j + 1) / (j + 1);
		}
	}

	double s{ x_term*B_term + A_term };
	return s;
}