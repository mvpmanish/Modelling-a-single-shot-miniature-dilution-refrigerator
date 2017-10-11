//cpp file to define functions declared in "vapour_pressure.hpp"

#include "vapour_pressure.hpp"

//------------------------------------------------------FUNCTION DEFINITIONS-----------------------------------------------------------------
//Function to calculate vapour pressure in [Pa]
double vap_pressure_He3(double T)
{

	double vap{ 0 };
	double temp{ 0 };
	for (int i{ 1 }; i < 6; i++)
	{
		temp += vp_He3_coeffs[i] * pow(T, i);
	}
	vap = vp_He3_coeffs[0] + temp;
	
	return vap;  //multiply by 100 to get pressure in [Pa] from [mBar]
}

//Function to calculate vapour pressure of Helium-4 in [Pa] from [Boreghi 1998]
double vap_pressure_He4(double T)
{
	double P{ exp(12.2440 - (59.83 / (8.31*T)) + (5 / 2)*log(T)) };
	return P;
}

//Function to write vapour pressure vs Temperature to a file
void write_vap_pressure_to_file(string filename, double start_temp, double end_temp, double temp_interval)
{
	ofstream(myfile); myfile.open(filename);
	myfile << "#T (Kelvin)  Vapour Pressure (Pa)" << endl;
	double d{ round((end_temp - start_temp)/temp_interval) };  //Take the upper integer for the number of intervals
	(int)( d ); //cast d to int
	for (int j{ 0 }; j < d + 1; j++)
	{
		double T{ start_temp + j*temp_interval };
		myfile << T << " " << vap_pressure_He3(T) << endl;
	}

	myfile.close(); cout << endl << "Succesfully read vapour pressure to " << filename << endl;
}