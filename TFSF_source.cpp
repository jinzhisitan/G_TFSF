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

	fout << setw(4) << nt*dt*1.0E6 << setw(14) << Ein(Isource) << setw(14) << Ein(-10) 
		<< setw(14) <<source(nt*dt-32.0/40.0/ c0)
		<< endl;
}

double TFSF::source(double time, double phase)
{
	////高空电磁脉冲源  IEC
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
	//升余弦函数
	double Ut = 1.0;
	double Tperiod = 2.0 * pi / omega;
	double t0 = Tperiod / 2.0;   //可以修改调整。
	if (time >= t0){
		Ut = 1.0;
	}
	else if (time >= 0.0){
		Ut = 0.5*(1.0 - cos(pi*time / t0)) ;
	}
	else {
		Ut = 0.0;
	}
	double fsource;
	fsource = Ut*sin(omega *time+phase);
	return fsource;


	//double f = 300.E6;
	//double tau = 2.0 / 150.0E6;
	//double t0 = 0.8*tau;
	//double fsource;
	//if (time >= 0) {
	//	fsource = sin(2.0 * pi *f* time)  *  exp(-4.0 *pi * (time - t0)* (time - t0) / tau/tau);
	//}
	//else {
	//	fsource = 0.0;
	//}
	//return fsource;
}