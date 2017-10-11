//Manish Patel and Adriano Parisi
//Last date modified 25/10/2016
//This is a program to model a single-shot miniature dilution refrigerator with time steps for the process of cooling

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>

#include "direct.h"
#include "mixing_chamber.hpp"
#include "still.hpp"
#include "molar_volumes.hpp"
#include "heat_leaks.hpp"

//-------------------------------------------------------------CONSTANTS-----------------------------------------------------------

//Definition of constants
double x_t{ 0.674 }, T_t{ 0.867 };  //Tricritical concentration and temperature[K]
double h_at_x_t{ -0.01986 }, dh_dx{ 0.0609 };  //Constants for specific heat in two-phase region [J/mol]
double He3_mass{ 5.008*pow(10, -27) }, He4_mass{ 6.646*pow(10, -27) };  //He-3 and He-4 atomic mass [kg]
double N_A{ 6.0221409*pow(10, 23) };  //Avagadro's number [/mol]
double He4_density{ 145 };  //density of He4 in [kg/m^3]
double He3_density{ 59 };  //density of He3 in [kg/m^3]
double gamma_e{ 0.69e-3 };  //Units of [(J/mol)*K^-2], coefficient for electronic heat capacity of copper [A.Tari 2003]
double pi{ 3.141592654 };

//Vectors of constants for fit determining specific heat of helium-3/helium-4 in two phase region [Chaudhry 2009]
vector<double> A2_1{ 17.221625, 12.661883, -963.11275, 7532.3747, -27191.567, 44976.189 };  //T < 0.15K
vector<double> A2_2{ 21.395625, -112.07332, 466.88150, -867.12111, 748.45512, -253.38487 };  // T > 0.15K
vector<double> B2_1{ 17.074164, -13.057881, 74.625853, -7421.4648, 56440.812, -116022.79 };  //T < 0.15K
vector<double> B2_2{ 20.841016, -117.76463, 500.33736, -998.13371, 927.34209, -328.83486 };  //T >0.15K

//Vector for vapour pressure curve constants for He-3
vector <double> vp_He3_coeffs{ -1.98e2,	1.9e3, -6.7e3, 1.07e4, -7.96e3, 3.44e3 }; // units of [Pa]

//Heat leak constants
double sb_const{ 5.6704e-8 };
double e_cu{ 0.7 }, e_st{ 0.5 };
double cf_leg_area{ 3.02e-6 };  //Area of one carbon-fibre support leg [m^2]
double p1_cf{ 8.8564e-2 }, p2_cf{ -1.7132e-6 };  //Constants which are from linear fit of heat leak vs temperature difference p1 in [Wm^-2K^-1] and p2 in [Wm^-2]

double still_temp{ 0.6 };  //Temperature of the still [K]

int main()
{
	ofstream datafile2; datafile2.open("still_temp_vs_cool_pwr_at_100mK_2uw_load.txt");
		cout << "Sup";
		//for (int g{ 0 }; g < 101; g++)
		//{
			//Initial paramters
			double T{ 0.1 };  //Starting temperatue of mixing chamber [Kelvin]
			double mc_pressure{ vap_pressure_He3(T) };  //Starting pressure in mixing chamber [Pa]
			double n3_dot{ 0 };  //Molar flow rate [mol/s], constant for now
			double n4_dot{ 0 };  //Molar flow rate of Helium-4 into mixing chamber due to pressure from Helium-3 leaving chamber [mol/s]
			double x{ 0.5 };  //Initial concentrations of helium-3 in mixture
			double q_dot;  //Cooling power [Watts]
			double mc_volume{ 8.67e-6 };  //Volume of mixing chamber [m^3]
			double molar_volume{ (1 - x)*He4_molar_vol(T, mc_pressure) + x*He3_molar_vol(T, mc_pressure) + x*(1 - x)*molar_vol_correction(T, mc_pressure, x) };  //Molar volume of the two-phase mixture of He-3 and He-4
			double n{ mc_volume / molar_volume };  //Amount of mixture at the start [mols]
			double n3{ n*x };  //Amount of helium-3 in system
			double q_load{ 2e-6};  //Heat load supplying a continuous heat supply to system in [W]
			double n_cu{ 4.474 };  //Number of mols of copper that makes up mixing chamber and bottom plate

			double still_volume{ 5.61e-5 };  //Volume of the still [m^3]
			//still_temp = 0.50+ 0.002*g;
			double start_hold_time{ 0 }, end_hold_time{ 0 }, counter{ 0 }, counter1{ 0 };  //Start and end times that mixing chamber is below T = 100mK

			mixing_chamber mc(T, mc_volume, n3, n - n3);  //Make mixing chamber for initial conditions
			still st(still_temp, still_volume);

			//Datafile name using strings
			ostringstream strs;
			strs << still_temp;
			ofstream datafile;
			string datafile_name("output_data_" + strs.str() + ".txt");
			//datafile.open(datafile_name);

			write_vap_pressure_to_file("new_vap_press_data.txt", 0.01, 0.9, 0.01);
			double d{ round((0.9 - 0.01) / 0.01) };  //Take the upper integer for the number of intervals
			(int)(d); //cast d to int
			ofstream datafilemolflow; datafilemolflow.open("new_mol_flow_vs_T.txt");
			datafilemolflow << "#T (Kelvin)  Molar flow rate (mol/s)" << endl;
			for (int i{ 0 }; i < d + 1; i++)
			{
				double Ti{ 0.01 + 0.01*i };
				datafilemolflow << Ti << " " << st.molecular_flow(Ti, vap_pressure_He3(Ti), 3.016e-3, 6.05e-3, 175e-3) << endl;
			}
			datafilemolflow.close();

			/*//------------------------------Steady-state: still_temp vs q_dot----------------------------------
			ofstream datafile3; datafile3.open("still_temp_vs_qdot_andy_correction.txt");
			for (int k{ 0 }; k < 326; k++)
			{
				still_temp = 0.475;
				//Calculate q_dot
				cout << k << "  ";
				still_temp = still_temp + 0.001*k;
				cout << still_temp << " ";
				still st(still_temp, still_volume);  
				q_dot = st.get_n3_dot()*84*pow(T,3);

				//Read in data
				if (counter == 0)
				{
					datafile3 << "%T_still [K]" << "  Cooling power q_dot [W]" << endl;
					counter++;
				}
				cout << still_temp << " " << q_dot << endl;
				datafile3 << still_temp << " " << q_dot << endl;
				
			} 
			datafile3.close()*/
			
			//---------------------Steady-state: heat load vs hold time-------------------------
			ofstream datafile4; datafile4.open("heat_load_vs_hold_time_vs_T_still.txt");
			for (int i{ 0 }; i < 51; i++)
			{
				if (i >= 0 && i <= 10)
				{
					q_dot = 1 + 0.15*i;
				}
				else
				{
					i = i - 10;
					q_dot = 2.5 + 0.5*i;
					i += 10;
				}
				
				still_temp = -4.3e-12*pow(q_dot, 6) + 1.3e-9*pow(q_dot, 5) - 1.7e-7*pow(q_dot, 4) + 1.2e-5*pow(q_dot, 3) - 0.00049*pow(q_dot, 2) + 0.016*q_dot + 0.43;
				still st1(still_temp, still_volume);
				n3_dot = st1.get_n3_dot();
				double t_hold{ n3 / n3_dot };
				cout << "q_dot=" << q_dot << " => " << "T_still=" << still_temp << " => t_hold=" << t_hold << endl;

				//Read in data
				if (counter == 0)
				{
					datafile4 << "%Heat load [W]" << "  Hold time [s]" << "T_still [K]" << endl;
					counter++;
				}

				datafile4 << q_dot << " " << t_hold << " " << still_temp << endl;
			} 
			datafile4.close();

			
			//------------------------------------Modelling cool down-------------------------------------
			/*//Put in for loop with t as time steps [seconds]
			for (int t1 = 0; t1 < 100000; t1++)
			{
				
				//Cooling power (i.e. number of Joules taken out of system in one second)
				q_dot = st.get_n3_dot()*mc.delta_h() - q_load; //- q_leak(T);

				//Read data into output_data
				if (datafile.is_open())
				{
					datafile << t1 << "  " << q_dot << "  " << T << "  " << mc.mol_He3() << "  " << n3_dot << "\n";
				}

				//Decrease energy of the mixing chmaber by taking out cooling power number of joules in 1s
				mc.lower_energy(q_dot);
				cout << endl << t1 << "s => T = " << mc.temp() << "K n3 =" << mc.mol_He3() << " n4 = " << mc.mol_He4() << " Qdot =" << q_dot << "W";
				//cout << endl << "C = " << mc.specific_heat();

				//Decrease amount of Helium-3 in mixing chamber
				n3_dot = st.get_n3_dot();
				mc.lower_n3(n3_dot);
				//Increase amount of Helium-4 in mixing chamber
				mc_pressure = mc.pressure();
				T = mc.temp();  //T has changed
				n4_dot = n3_dot*He3_molar_vol(T, mc_pressure) / He4_molar_vol(T, mc_pressure);
				mc.increase_n4(n4_dot);

				//Set the start of the hold time
				if (mc.temp() <= 100e-3 && counter == 0)  //counter ensures that in the loop this statement is accessed only once
				{
					start_hold_time = t1;
					datafile2 << still_temp << " " << q_dot << endl;
					counter++;
				}

				if (mc.mol_He3() <= 0)
				{
					cout << endl << "Helium-3 has run out!" << endl;
					double min_T{ T };

					//Helium-3 has run out and let the system warm up for t2 seconds
					for (int t2 = 0; t2 < 3600; t2++)
					{
						double q{ -q_load }; //- q_leak(T) };  //Heat is flowing into system now
						mc.lower_energy(q);
						T = mc.temp();
						cout << endl << t1 + t2 << "s => T = " << mc.temp() << "K";

						//Read data to datafile
						if (datafile.is_open())
						{
							datafile << t1 + t2 << " " << T << "\n";
						} 

						if (T >= 100e-3 && counter == 1)  //Set the end of the hold time
						{
							end_hold_time = t1 + t2;
							counter++;
						}
					}

					double hold_time{ (end_hold_time - start_hold_time) / 60 };  //Hold time in minutes

					cout << endl << "Minimum temperature of mixing chamber = " << min_T << " K" << endl
						<< "Minimum cooling power = " << q_dot << " Watts" << endl
						<< "Cooling time = " << t1 / 60 << " mins" << endl
						<< "Hold time = " << hold_time << " mins" << endl;
					
					datafile.close();
					break;
				}
			}
			
		//} */
	datafile2.close();
	cout << endl << "Finished. (Press enter to exit)";
	cin.ignore(); //User has to press enter to exit at end of program
}