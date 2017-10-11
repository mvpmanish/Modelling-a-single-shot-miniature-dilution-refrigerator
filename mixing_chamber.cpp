//cpp file to define the functions declared in "mixing_chamber.hpp"

#include "mixing_chamber.hpp"

//--------------------------------------------TWO-PHASE FUNCTIONS---------------------------------------------------
//Dilute phase He-3 concentration
double two_phase_x_dilute(double &T)
{
	double x_dil{ 0 };
	if (T < 0.2)
	{
		x_dil = 0.066 - 0.027505*T + 0.68694*pow(T, 2);
		return x_dil;
	}
	x_dil = (-0.209148*(T - T_t) / (T - T_t - 0.080280)) + (0.960222*(T - T_t)) + (0.549920*pow(T - T_t, 2));
	return x_dil;
}

//Concentrated phase He-3 concentration
double two_phase_x_conc(double &T)
{
	double x_conc{ 0 };
	if (T < 0.2)
	{
		x_conc = 1;
		return x_conc;
	}
	x_conc = (0.316170*pow(T - T_t, 3)) - (0.180743*pow(T - T_t, 2)) - (0.746805*(T - T_t));
	return x_conc;
}

//Specific heat in the two phase region [J/(K*kg)]
double two_phase_specific_heat(double &T)
{
	double C{ 0 };  //Calculate specific heat by adding terms of T^i using vectors of constants
	if (T < 0.15)  //T < 0.15K
	{
		for (int i = 0; i < 6; i++)
		{
			C += A2_1[i] * pow(T, i);
		}
	}
	else  //T >= 0.15K
	{
		for (int i = 0; i < 6; i++)
		{
			C += A2_2[i] * pow(T, i);
		}
	}

	return C;
}

//--------------------------------------------MIXING CAHMBER FUNCTION DEFINITIONS----------------------------------------------

//Parameterised constrcutor
mixing_chamber::mixing_chamber(double Ti, double Vi, double n3i, double n4i)
{
	T = Ti; V = Vi; n3 = n3i; n4 = n4i;

	x_dil = two_phase_x_dilute(Ti); //dilute phase concentration of He-3
	x_conc = two_phase_x_conc(Ti);  //concentrated phase concentration of He-3

	C = two_phase_specific_heat(Ti);

	P = vap_pressure_He3(Ti);  //Vapour pressure of Helium-3 (don't need to include He-4 as it is a factor of 10^3 times smaller at these temperatures)
}

double mixing_chamber::h_dil()  //Enthalpy of dilute phase
{
	double h2{ 0 };
	if (T < 0.15)  //Use constants A2_1[i] and B2_1[i]
	{
		for (int i = 0; i < 6; i++)
		{
			h2 += (B2_1[i] * (x_dil - x_t) + A2_1[i]) * pow(T, i + 2) / (i + 2);
		}

		h2 += dh_dx + h_at_x_t;
	}
	else  //T >= 0.15K, use constants A2_2[i] and B2_2[i]
	{
		for (int i = 0; i < 6; i++)
		{
			h2 += (B2_2[i] * (x_dil - x_t) + A2_2[i]) * pow(T, i + 2) / (i + 2);
		}

		h2 += dh_dx + h_at_x_t;
	}

	return h2;
}

double mixing_chamber::h_conc()
{
	double h2{ 0 };
	if (T < 0.15)  //Use constants A2_1[i] and B2_1[i]
	{
		for (int i = 0; i < 6; i++)
		{
			h2 += (B2_1[i] * (x_conc - x_t) + A2_1[i]) * pow(T, i + 2) / (i + 2);
		}

		h2 += dh_dx + h_at_x_t;
	}
	else  //T >= 0.15K, use constants A2_2[i] and B2_2[i]
	{
		for (int i = 0; i < 6; i++)
		{
			h2 += (B2_2[i] * (x_conc - x_t) + A2_2[i]) * pow(T, i + 2) / (i + 2);
		}

		h2 += dh_dx + h_at_x_t;
	}

	return h2;
}



