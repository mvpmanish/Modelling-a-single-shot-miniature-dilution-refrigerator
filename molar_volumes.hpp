//Header file to define the functions for the molar volume of pure He-4 and pure He-3 at a given temperature [K] and pressure [Pa] from Chaudhry 2009

#ifndef MOLAR_VOLUMES_HPP
#define MOLAR_VOLUMES_HPP

#include<vector>
#include<cmath>

using namespace std;

//---------------------------------------CONSTANTS------------------------------------
//Constants for calculating molar volume of He-4
vector<double> v40{ 2.757930e-5, -3.361585e-12, 1.602419e-18, -1.072604e-24, 7.979064e-31, -5.356076e-37, 2.703689e-43, -9.004790e-50, 1.725962e-56, -1.429411e-63 };  //Units of [(m^3/mol)*Pa^i]
vector<double> Delta{ 8.618, -7.487e-7, 5.308e-14 };  //Units of [K/(Pa^i)]
double A1{ 1.863e8 };  //Units [Pa^2/K^4]
double A2{ 7.664e6 };  //Units [Pa^(5/3)K^4]
double p_c{ -1.00170e6 };  //Units [Pa]
double B0{ 2.191 }, B1{ -1.702e-7 };  //B1 has units [Pa^-1]

//Constants for calculating molar volume of He-3
vector<vector<double>> v31{ { 40.723012, -1.3151948e-5, 0.0409498e-10 }, { -0.6614794, 0.1275125e-5, -0.0091931e-10 }, { 0.5542147, -0.1527959e-5, 0.0113764e-10 }, { 0.1430724, -0.0034712e-5, -0.0008515e-10 } }; //Matrix v31[i][j] in units [(cm^3)*(mol^-1)*(K^-i)(Pa^-j)]
double v320{ -0.2603492 }, v322{ -0.0051946e-10 }; //v320 units [mol/cm^3] and v322 units [mol*(cm^-3)*(Pa)]

//Constants for claculating correction term: use units of cm^3/mol and pressures in bar (10^5 Pa)
vector<vector<double>> vr10{ { 0.1389646, -2.5679051, 0.5113374, 3.9851834 }, { -4.933308, 40.217899, -83.295844, 16.616527 }, { -4.4736292, 21.103275, 42.535669, 32.247635 }, { -17.763600, -17.585800, -2.7419703, -21.867304 } };
vector<vector<double>> vr11{ { -0.0706271, 0.8984671, -0.1993514, 2.3961147 }, { 2.0430264, -13.570937, 18.303655, 2.3961147 }, { 1.9707327, 1.2376227, -9.101332, -18.726933 }, { 1.2116373, 0.3340713, 7.5264822, 7.8044314 } };
vector<vector<double>> vr12{ { 0.010287, -0.0599619, -0.0262424, 0.1356818 }, { -0.1641865, 0.972759, -0.9242861, -0.5426784 }, { -0.2765908, 0.2888217, -0.6437575, 2.2866511 }, { 0.0781768, -0.3641635, 0.5378945, -1.3449638 } };
vector<vector<double>> vr2{ { -1.1647447, -4.5046372, 10.15904, -11.353903 }, { -8.2289003, 24.648775, 31.640583, 4.2288385 }, { -5.7908683, -104.81653, 29.677981, -81.740459 }, { 62.107622, 8.1406542, 4.3197622, 32.623983 } };
double p_vr{ 1.7867623 };

//--------------------------------------------------FUNCTIONS--------------------------------------------------
//Molar volume of pure He-4 in [m^3/mol]
double He4_molar_vol(double &T, double &p)
{
	double v4_1{ 0 };
	for (int i{ 0 }; i < 10; i++)
	{
		v4_1 += v40[i] * pow(p, i);
	}

	double v4_2{ 0 };
	for (int i{ 0 }; i < 3; i++)
	{
		v4_2 += Delta[i] * pow(p, i); 
	}
	
	double v4_3{ 1 + ((A1 / pow(p - p_c, 2)) + (A2 / pow(p - p_c, 5 / 3)))*(pow(T, 4) / 4) };
	double v4_4{ v4_1*(v4_3 - (B0 + B1*p)*erfc(pow(v4_2 / T, 0.5))) };

	return v4_4;
}

//Molar volume of He-3 in [m^3/mol]
double He3_molar_vol(double &T, double &p)
{
	double v3_1{ 0 };
	for (int i{ 0 }; i < 4; i++)
	{
		for (int j{ 0 }; j < 3; j++)
		{
			v3_1 += v31[i][j] * pow(T, i) * pow(p, j);
		}
	}
	double v3_2{ v3_1 + (1 / (v322*pow(p, 2) + v320)) };

	return v3_2/1e6;  //divide by 10^6 to get molar volume in [m^3/mol] instead of [cm^3/mol]
}

//Correction term for molar volume in [m^3/mol]
double molar_vol_correction(double &T, double &p_i, double &x)
{
	double vr{ 0 };
	double p{ p_i/1e5 };  //Convert pressure into bar from Pa

	for (int i{ 0 }; i < 4; i++)
	{
		for (int j{ 0 }; j < 4; j++)
		{
			vr += pow(x, i)*pow(T, j)*(vr10[i][j] + vr11[i][j] * p + vr12[i][j] * p + vr2[i][j] / (p + p_vr));
		}
	}

	return vr/1e6;  //divide by 10^6 to get molar volume in [m^3/mol] instead of [cm^3/mol]
}

#endif