#include<complex>
#include"TFSF.h"

double TFSF::source(double time)
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

	////////IEC 0-30MHz
	//const double pi = acos(-1.0);
	//const complex<double>ci(0.0, 1.0);
	//complex<double>ftemp;
	//double Ttotal = 0.5e-6;
	//double Fresolution = 1.0 / Ttotal;
	//if (time >=0.0) {   //&&time<= Ttotal
	//	fsource = E0k*(1.0 / alpha - 1.0 / beta) / Ttotal;
	//	for (int k = 1; k <= 15; k++) {
	//		ftemp = E0k*(1.0 / (alpha + k * 2.0 * pi*Fresolution*ci)
	//			- 1.0 / (beta + k* 2.0 * pi*Fresolution*ci));
	//		fsource = fsource + 2.0 * real(ftemp*exp(k * 2.0 * pi*Fresolution*time*ci)) / Ttotal;
	//	}
	//}
	//else {
	//	fsource = 0.0;
	//}
	//return fsource;

	////////高斯
	//double fsource;
	//if (time >= 0) {
	//	fsource = sin(5.34E9 *time)*exp(-4.0*pi*
	//		(time - 3.53E-9) / 3.13E-9*(time - 3.53E-9) / 3.13E-9);
	//}
	//else {
	//	fsource = 0.0;
	//}
	//return fsource;

	/////单一频率
	double omega = 2.0 * pi*300.0E6;
	double fsource;
	if (time >= 0) {
		fsource = sin(omega *time);
	}
	else {
		fsource = 0.0;
	}
	return fsource;
}