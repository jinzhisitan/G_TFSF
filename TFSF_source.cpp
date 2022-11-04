#include<cmath>
#include<iostream>
#include <fstream>
#include<iomanip>
#include"physical_constants.h"
#include"TFSF.h"

using namespace std;

std::string fname = "zz_TFSF.dat";
ofstream fout(fname);

void TFSF::addSource(int nt)
{
	Ein(Isource) = source(nt*dt);

	fout << scientific;
	fout.precision(4);
	fout << setw(4) << nt << setw(14) << Ein(Isource) << setw(14) << Ein(0) << endl;
}

double TFSF::source(double time)
{
	////¸ß¿Õµç´ÅÂö³åÔ´  IEC
	//int option = 3;
	//double E0k, alpha, beta, fsource;
	//// Publication in 1976
	//if (option == 1) {
	//	E0k = 1.04*5.0E4;
	//	alpha = 1.5E6;
	//	beta = 2.6E8;
	//}
	//else if (option == 2) {//Bell Lab.
	//	E0k = 1.05*5.0E4;
	//	alpha = 4.0E6;
	//	beta = 4.76E8;
	//}
	//else { //  IEC
	//	E0k = 1.3*5.0E4;
	//	alpha = 40.0E6;
	//	beta = 6.0E8;
	//}
	//if (time >= 0) {
	//	fsource = E0k * (exp(-alpha * time) - exp(-beta * time));
	//}
	//else {
	//	fsource = 0.0;
	//}
	//return fsource;

	double omega = 2.0 * pi*300.E6;
	double fsource;
	if (time >= 0) {
		fsource = sin(omega *time);
	}
	else {
		fsource = 0.0;
	}
	return fsource;


	//double f = 300.E6;
	//double tau = 2.0 / 100.0E6;
	//double t0 = 0.6*tau;
	//double fsource;
	//if (time >= 0) {
	//	fsource = sin(2.0 * pi *f* time)  *  exp(-4.0 *pi * (time - t0)* (time - t0) / tau/tau);
	//}
	//else {
	//	fsource = 0.0;
	//}
	//return fsource;
}