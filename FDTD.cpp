#include<string>
#include <fstream>
#include"FDTD.h"

FDTD::FDTD()
{
	std::cout << "class FDTD: void constructor\n";
}

FDTD::FDTD(int tmax_, int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_, int xPML1, int xPML2,
	int yPML1, int yPML2, int zPML1, int zPML2, double dt_, double dx_, double dy_, double dz_, double epsr_[], double sig_[],
	double alpha, double thi, double phi, int IncStart, int IncEnd, int ItMin, int ItMax, int JtMin, int JtMax, int KtMin, int KtMax)
	:Einc(dx_, dy_, dz_, alpha, thi, phi,  dt_, IncStart, IncEnd, ItMin, ItMax, JtMin, JtMax, KtMin, KtMax, IsMin_, IsMax_, JsMin_, JsMax_, KsMin_, KsMax_)
{
	IsMin = IsMin_; IsMax = IsMax_;
	JsMin = JsMin_; JsMax = JsMax_;
	KsMin = KsMin_; KsMax = KsMax_;

	nxPML_1 = xPML1; nxPML_2 = xPML2;
	nyPML_1 = yPML1; nyPML_2 = yPML2;
	nzPML_1 = zPML1; nzPML_2 = zPML2;
	tmax = tmax_;

	dx = dx_;
	dy = dy_;
	dz = dz_;
	dt = dt_;

	for (int i = 0; i < MediaNo; i++) {
		epsr[i] = epsr_[i];
		sig[i] = sig_[i];
	}

	Imin = IsMin_ - xPML1;
	Imax = IsMax_ + xPML2;
	Jmin = JsMin_ - yPML1;
	Jmax = JsMax_ + yPML2;
	Kmin = KsMin_ - zPML1;
	Kmax = KsMax_ + zPML2;

	Ex.Allocate(Imin, Imax - 1, Jmin, Jmax, Kmin, Kmax);
	Ey.Allocate(Imin, Imax, Jmin, Jmax - 1, Kmin, Kmax);
	Ez.Allocate(Imin, Imax, Jmin, Jmax, Kmin, Kmax - 1);
	Hx.Allocate(Imin, Imax, Jmin, Jmax - 1, Kmin, Kmax - 1);
	Hy.Allocate(Imin, Imax - 1, Jmin, Jmax, Kmin, Kmax - 1);
	Hz.Allocate(Imin, Imax - 1, Jmin, Jmax - 1, Kmin, Kmax);

	ob.Allocate(Imin, Imax - 1, Jmin, Jmax - 1, Kmin, Kmax - 1);

	den_ex.Allocate(Imin, Imax), den_ey.Allocate(Jmin, Jmax), den_ez.Allocate(Kmin, Kmax);
	den_hx.Allocate(Imin, Imax - 1), den_hy.Allocate(Jmin, Jmax - 1), den_hz.Allocate(Kmin, Kmax - 1);

	/////////////////////////////////////////////////////////////////////////////////
	///////////////////////CPML
	/////////////////////////////////////////////////////////////////////////////////
	psi_Exy_1.Allocate(Imin, Imax - 1, Jmin, JsMin, Kmin, Kmax); psi_Exy_2.Allocate(Imin, Imax - 1, JsMax, Jmax, Kmin, Kmax);
	psi_Exz_1.Allocate(Imin, Imax - 1, Jmin, Jmax, Kmin, KsMin), psi_Exz_2.Allocate(Imin, Imax - 1, Jmin, Jmax, KsMax, Kmax);
	psi_Eyx_1.Allocate(Imin, IsMin, Jmin, Jmax - 1, Kmin, Kmax), psi_Eyx_2.Allocate(IsMax, Imax, Jmin, Jmax - 1, Kmin, Kmax);
	psi_Eyz_1.Allocate(Imin, Imax, Jmin, Jmax - 1, Kmin, KsMin), psi_Eyz_2.Allocate(Imin, Imax, Jmin, Jmax - 1, KsMax, Kmax);
	psi_Ezx_1.Allocate(Imin, IsMin, Jmin, Jmax, Kmin, Kmax - 1), psi_Ezx_2.Allocate(IsMax, Imax, Jmin, Jmax, Kmin, Kmax - 1);
	psi_Ezy_1.Allocate(Imin, Imax, Jmin, JsMin, Kmin, Kmax - 1), psi_Ezy_2.Allocate(Imin, Imax, JsMax, Jmax, Kmin, Kmax - 1);

	psi_Hxy_1.Allocate(Imin, Imax, Jmin, JsMin - 1, Kmin, Kmax - 1), psi_Hxy_2.Allocate(Imin, Imax, JsMax, Jmax - 1, Kmin, Kmax - 1);
	psi_Hxz_1.Allocate(Imin, Imax, Jmin, Jmax - 1, Kmin, KsMin - 1), psi_Hxz_2.Allocate(Imin, Imax, Jmin, Jmax - 1, KsMax, Kmax - 1);
	psi_Hyx_1.Allocate(Imin, IsMin - 1, Jmin, Jmax, Kmin, Kmax - 1), psi_Hyx_2.Allocate(IsMax, Imax - 1, Jmin, Jmax, Kmin, Kmax - 1);
	psi_Hyz_1.Allocate(Imin, Imax - 1, Jmin, Jmax, Kmin, KsMin - 1), psi_Hyz_2.Allocate(Imin, Imax - 1, Jmin, Jmax, KsMax, Kmax - 1);
	psi_Hzx_1.Allocate(Imin, IsMin - 1, Jmin, Jmax - 1, Kmin, Kmax), psi_Hzx_2.Allocate(IsMax, Imax - 1, Jmin, Jmax - 1, Kmin, Kmax);
	psi_Hzy_1.Allocate(Imin, Imax - 1, Jmin, JsMin - 1, Kmin, Kmax), psi_Hzy_2.Allocate(Imin, Imax - 1, JsMax, Jmax - 1, Kmin, Kmax);

	be_x_1.Allocate(Imin, IsMin), ce_x_1.Allocate(Imin, IsMin);
	be_x_2.Allocate(IsMax, Imax), ce_x_2.Allocate(IsMax, Imax);
	be_y_1.Allocate(Jmin, JsMin), ce_y_1.Allocate(Jmin, JsMin);
	be_y_2.Allocate(JsMax, Jmax), ce_y_2.Allocate(JsMax, Jmax);
	be_z_1.Allocate(Kmin, KsMin), ce_z_1.Allocate(Kmin, KsMin);
	be_z_2.Allocate(KsMax, Kmax), ce_z_2.Allocate(KsMax, Kmax);

	bh_x_1.Allocate(Imin, IsMin - 1), ch_x_1.Allocate(Imin, IsMin - 1);
	bh_x_2.Allocate(IsMax, Imax - 1), ch_x_2.Allocate(IsMax, Imax - 1);
	bh_y_1.Allocate(Jmin, JsMin - 1), ch_y_1.Allocate(Jmin, JsMin - 1);
	bh_y_2.Allocate(JsMax, Jmax - 1), ch_y_2.Allocate(JsMax, Jmax - 1);
	bh_z_1.Allocate(Kmin, KsMin - 1), ch_z_1.Allocate(Kmin, KsMin - 1);
	bh_z_2.Allocate(KsMax, Kmax - 1), ch_z_2.Allocate(KsMax, Kmax - 1);
}


void FDTD::initCoeficients()
{
	//计算cell介电常数
	IsMedia();

///////////////////////////////////////////////////////////////////////////
	double temp;
	//FILL IN UPDATING COEFFICIENTS
	DA = 1.0;
	DB = dt / mu0;
	double epsr_eff, sig_eff;
	for (int i = 0; i < MediaNo; i++) {
		for (int j = 0; j < MediaNo; j++) {
			for (int k = 0; k < MediaNo; k++) {
				for (int l = 0; l < MediaNo; l++) {
					epsr_eff = 0.25*(epsr[i] + epsr[j] + epsr[k] + epsr[l]);    //4个网格介质的平均
					sig_eff = 0.25*(sig[i] + sig[j] + sig[k] + sig[l]);
					temp = sig_eff*dt / (2.0*eps0*epsr_eff);
					CA[i][j][k][l] = (1 - temp) / (1 + temp);      //Eq(2-2-17)
					CB[i][j][k][l] = dt / eps0 / (epsr_eff + sig_eff*dt / (2.0*eps0)); //Eq(2-2-18)
				}
			}	
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//FILL IN DENOMINATORS FOR FIELD UPDATES
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (int i = Imin; i <= Imax-1; i++) {
		den_hx(i) = 1.0 / dx;
	}
	for (int j = Jmin; j <= Jmax-1; j++) {
		den_hy(j) = 1.0 / dy;
	}
	for (int k = Kmin; k <= Kmax-1; k++) {
		den_hz(k) = 1.0 / dz;
	}
	for (int i = Imin; i <= Imax; i++) {
		den_ex(i) = 1.0 / dx;
	}
	for (int j = Jmin; j <= Jmax; j++) {
		den_ey(j) = 1.0 / dy;
	}
	for (int k = Kmin; k <= Kmax; k++) {
		den_ez(k) = 1.0 / dz;
	}


//////////////////////////////////////////////////////////////////////////////////
	//Specify the CPML Order and Other Parameters
	int m = 4, ma = 1;
	//[1]Eq(5-4-50), [2]Eq(7.66),Eq(7.67)
	double sig_x_max = (m + 1) / (dx * 150 * pi*sqrt(epsr[0]));
	double sig_y_max = sig_x_max;
	double sig_z_max = sig_x_max;
	//[1]P107  alpha=0~0.05 
	double alpha_x_max = 0.0, alpha_y_max = alpha_x_max, alpha_z_max = alpha_x_max;
	//[1]P107
	double kappa_x_max = 1.0, kappa_y_max = kappa_x_max, kappa_z_max = kappa_x_max;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//SET CPML PARAMETERS IN EACH DIRECTION
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int i, j, k;
	//  Ex  PML1
	double sige_x_PML_1, alphae_x_PML_1, kappae_x_PML_1;
	for (i = Imin; i <= IsMin; i++) {
		//[1]Eq(5-4-49)   [2]Eq(7.60a)
		sige_x_PML_1 = sig_x_max * pow((IsMin - i) / double(nxPML_1), m);
		//[1]Eq(5-4-52)   [2]Eq(7.79)  
		alphae_x_PML_1 = alpha_x_max*pow((i - Imin) / double(nxPML_1), ma);     //修改
																				//[1]Eq(5-4-51)   [2]Eq(7.60b)
		kappae_x_PML_1 = 1.0 + (kappa_x_max - 1.0)*pow((IsMin - i) / double(nxPML_1), m);

		be_x_1(i) = exp(-(sige_x_PML_1 / kappae_x_PML_1 + alphae_x_PML_1)  *dt / eps0);
		if ((sige_x_PML_1 == 0.0) && (alphae_x_PML_1 == 0.0) && (i == IsMin)) {
			ce_x_1(i) = 0.0;
		}
		else {
			ce_x_1(i) = sige_x_PML_1*(be_x_1(i) - 1.0) /
				(sige_x_PML_1 + kappae_x_PML_1*alphae_x_PML_1) / kappae_x_PML_1;
		}
		ce_x_1(i) = ce_x_1(i) / dx;
		//FILL IN DENOMINATORS FOR FIELD UPDATES
		den_ex(i) = 1.0 / (kappae_x_PML_1*dx);
	}

	//  Hx  PML1
	double sigh_x_PML_1, alphah_x_PML_1, kappah_x_PML_1;
	for (i = Imin; i <= IsMin - 1; i++) {
		sigh_x_PML_1 = sig_x_max * pow((IsMin - i - 0.5) / double(nxPML_1), m);
		alphah_x_PML_1 = alpha_x_max*pow((i + 0.5 - Imin) / double(nxPML_1), ma);
		kappah_x_PML_1 = 1.0 + (kappa_x_max - 1.0)* pow((IsMin - i - 0.5) / double(nxPML_1), m);

		bh_x_1(i) = exp(-(sigh_x_PML_1 / kappah_x_PML_1 + alphah_x_PML_1)*dt / eps0);
		ch_x_1(i) = sigh_x_PML_1*(bh_x_1(i) - 1.0) /
			(sigh_x_PML_1 + kappah_x_PML_1*alphah_x_PML_1) / kappah_x_PML_1;

		ch_x_1(i) = ch_x_1(i) / dx;
		den_hx(i) = 1.0 / (kappah_x_PML_1*dx);
	}


	//  Ex  PML2
	double sige_x_PML_2, alphae_x_PML_2, kappae_x_PML_2;
	for (i = IsMax; i <= Imax; i++)
	{
		sige_x_PML_2 = sig_x_max * pow((i - IsMax) / double(nxPML_2), m);
		alphae_x_PML_2 = alpha_x_max*pow((Imax - i) / double(nxPML_2), ma);
		kappae_x_PML_2 = 1.0 + (kappa_x_max - 1.0)*pow((i - IsMax) / double(nxPML_2), m);

		be_x_2(i) = exp(-(sige_x_PML_2 / kappae_x_PML_2 + alphae_x_PML_2)*dt / eps0);
		if ((sige_x_PML_2 == 0.0) && (alphae_x_PML_2 == 0.0) && (i == IsMax)) {
			ce_x_2(i) = 0.0;
		}
		else {
			ce_x_2(i) = sige_x_PML_2*(be_x_2(i) - 1.0) /
				(sige_x_PML_2 + kappae_x_PML_2*alphae_x_PML_2) / kappae_x_PML_2;
		}
		ce_x_2(i) = ce_x_2(i) / dx;
		den_ex(i) = 1.0 / (kappae_x_PML_2*dx);
	}

	//  Hx  PML2
	double sigh_x_PML_2, alphah_x_PML_2, kappah_x_PML_2;
	for (i = IsMax; i <= Imax - 1; i++)
	{
		sigh_x_PML_2 = sig_x_max * pow((i + 0.5 - IsMax) / double(nxPML_2), m);
		alphah_x_PML_2 = alpha_x_max*pow((Imax - i - 0.5) / double(nxPML_2), ma);
		kappah_x_PML_2 = 1.0 + (kappa_x_max - 1.0)*pow((i + 0.5 - IsMax) / double(nxPML_2), m);

		bh_x_2(i) = exp(-(sigh_x_PML_2 / kappah_x_PML_2 + alphah_x_PML_2)*dt / eps0);
		ch_x_2(i) = sigh_x_PML_2*(bh_x_2(i) - 1.0) /
			(sigh_x_PML_2 + kappah_x_PML_2*alphah_x_PML_2) / kappah_x_PML_2;

		ch_x_2(i) = ch_x_2(i) / dx;
		den_hx(i) = 1.0 / (kappah_x_PML_2*dx);
	}

	//  Ey  PML1  
	double sige_y_PML_1, alphae_y_PML_1, kappae_y_PML_1;
	for (j = Jmin; j <= JsMin; j++)
	{
		sige_y_PML_1 = sig_y_max * pow((JsMin - j) / double(nyPML_1), m);
		alphae_y_PML_1 = alpha_y_max*pow((j - Jmin) / double(nyPML_1), ma);
		kappae_y_PML_1 = 1.0 + (kappa_y_max - 1.0)*pow((JsMin - j) / double(nyPML_1), m);

		be_y_1(j) = exp(-(sige_y_PML_1 / kappae_y_PML_1 + alphae_y_PML_1)*dt / eps0);
		if ((sige_y_PML_1 == 0.0) && (alphae_y_PML_1 == 0.0) && (j == JsMin)) {
			ce_y_1(j) = 0.0;
		}
		else {
			ce_y_1(j) = sige_y_PML_1*(be_y_1(j) - 1.0) /
				(sige_y_PML_1 + kappae_y_PML_1*alphae_y_PML_1) / kappae_y_PML_1;
		}
		ce_y_1(j) = ce_y_1(j) / dy;
		den_ey(j) = 1.0 / (kappae_y_PML_1*dy);
	}

	//  Hy  PML1
	double sigh_y_PML_1, alphah_y_PML_1, kappah_y_PML_1;
	for (j = Jmin; j <= JsMin - 1; j++)
	{
		sigh_y_PML_1 = sig_y_max * pow((JsMin - j - 0.5) / double(nyPML_1), m);
		alphah_y_PML_1 = alpha_y_max*pow((j + 0.5 - Jmin) / double(nyPML_1), ma);
		kappah_y_PML_1 = 1.0 + (kappa_y_max - 1.0)*pow((JsMin - j - 0.5) / double(nyPML_1), m);

		bh_y_1(j) = exp(-(sigh_y_PML_1 / kappah_y_PML_1 + alphah_y_PML_1)*dt / eps0);
		ch_y_1(j) = sigh_y_PML_1*(bh_y_1(j) - 1.0) /
			(sigh_y_PML_1 + kappah_y_PML_1*alphah_y_PML_1) / kappah_y_PML_1;

		ch_y_1(j) = ch_y_1(j) / dy;
		den_hy(j) = 1.0 / (kappah_y_PML_1*dy);
	}

	//  Ey PML2
	double sige_y_PML_2, alphae_y_PML_2, kappae_y_PML_2;
	for (j = JsMax; j <= Jmax; j++)
	{
		sige_y_PML_2 = sig_y_max * pow((j - JsMax) / double(nyPML_2), m);
		alphae_y_PML_2 = alpha_y_max*pow((Jmax - j) / double(nyPML_2), ma);
		kappae_y_PML_2 = 1.0 + (kappa_y_max - 1.0)*pow((j - JsMax) / double(nyPML_2), m);

		be_y_2(j) = exp(-(sige_y_PML_2 / kappae_y_PML_2 + alphae_y_PML_2)*dt / eps0);
		if ((sige_y_PML_2 == 0.0) && (alphae_y_PML_2 == 0.0) && (j == JsMax)) {
			ce_y_2(j) = 0.0;
		}
		else {
			ce_y_2(j) = sige_y_PML_2*(be_y_2(j) - 1.0)
				/ (sige_y_PML_2 + kappae_y_PML_2*alphae_y_PML_2) / kappae_y_PML_2;
		}
		ce_y_2(j) = ce_y_2(j) / dy;
		den_ey(j) = 1.0 / (kappae_y_PML_2*dy);
	}
	//  Hy PML2
	double sigh_y_PML_2, alphah_y_PML_2, kappah_y_PML_2;
	for (j = JsMax; j <= Jmax - 1; j++)
	{
		sigh_y_PML_2 = sig_y_max * pow((j + 0.5 - JsMax) / double(nyPML_2), m);
		alphah_y_PML_2 = alpha_y_max*pow((Jmax - j - 0.5) / double(nyPML_2), m);
		kappah_y_PML_2 = 1.0 + (kappa_y_max - 1.0)*pow((j + 0.5 - JsMax) / double(nyPML_2), m);

		bh_y_2(j) = exp(-(sigh_y_PML_2 / kappah_y_PML_2 + alphah_y_PML_2)*dt / eps0);
		ch_y_2(j) = sigh_y_PML_2*(bh_y_2(j) - 1.0) /
			(sigh_y_PML_2 + kappah_y_PML_2*alphah_y_PML_2) / kappah_y_PML_2;

		ch_y_2(j) = ch_y_2(j) / dy;
		den_hy(j) = 1.0 / (kappah_y_PML_2*dy);
	}

	//  Ez PML1
	double sige_z_PML_1, alphae_z_PML_1, kappae_z_PML_1;
	for (k = Kmin; k <= KsMin; k++) {
		sige_z_PML_1 = sig_z_max * pow((KsMin - k) / double(nzPML_1), m);
		alphae_z_PML_1 = alpha_z_max*pow((k - Kmin) / double(nzPML_1), ma);
		kappae_z_PML_1 = 1.0 + (kappa_z_max - 1.0)*pow((KsMin - k) / double(nzPML_1), m);

		be_z_1(k) = exp(-(sige_z_PML_1 / kappae_z_PML_1 + alphae_z_PML_1)*dt / eps0);
		if ((sige_z_PML_1 == 0.0) && (alphae_z_PML_1 == 0.0) && (k == KsMin)) {
			ce_z_1(k) = 0.0;
		}
		else {
			ce_z_1(k) = sige_z_PML_1*(be_z_1(k) - 1.0) /
				(sige_z_PML_1 + kappae_z_PML_1*alphae_z_PML_1) / kappae_z_PML_1;
		}
		ce_z_1(k) = ce_z_1(k) / dz;
		den_ez(k) = 1.0 / (kappae_z_PML_1*dz);
	}
	//  Hz PML1
	double sigh_z_PML_1, alphah_z_PML_1, kappah_z_PML_1;
	for (k = Kmin; k <= KsMin - 1; k++)
	{
		sigh_z_PML_1 = sig_z_max * pow((KsMin - k - 0.5) / double(nzPML_1), m);
		alphah_z_PML_1 = alpha_z_max*pow((k + 0.5 - Kmin) / double(nzPML_1), ma);
		kappah_z_PML_1 = 1.0 + (kappa_z_max - 1.0)*pow((KsMin - k - 0.5) / double(nzPML_1), m);

		bh_z_1(k) = exp(-(sigh_z_PML_1 / kappah_z_PML_1 + alphah_z_PML_1)*dt / eps0);
		ch_z_1(k) = sigh_z_PML_1*(bh_z_1(k) - 1.0) /
			(sigh_z_PML_1 + kappah_z_PML_1*alphah_z_PML_1) / kappah_z_PML_1;

		ch_z_1(k) = ch_z_1(k) / dz;
		den_hz(k) = 1.0 / (kappah_z_PML_1*dz);
	}

	//  Ez PML2
	double sige_z_PML_2, alphae_z_PML_2, kappae_z_PML_2;
	for (k = KsMax; k <= Kmax; k++)
	{
		sige_z_PML_2 = sig_z_max * pow((k - KsMax) / double(nzPML_2), m);
		alphae_z_PML_2 = alpha_z_max*pow((Kmax - k) / double(nzPML_2), m);
		kappae_z_PML_2 = 1.0 + (kappa_z_max - 1.0)*pow((k - KsMax) / double(nzPML_2), m);

		be_z_2(k) = exp(-(sige_z_PML_2 / kappae_z_PML_2 + alphae_z_PML_2)*dt / eps0);
		if ((sige_z_PML_2 == 0.0) && (alphae_z_PML_2 == 0.0) && (k == KsMax)) {
			ce_z_2(k) = 0.0;
		}
		else {
			ce_z_2(k) = sige_z_PML_2*(be_z_2(k) - 1.0) /
				(sige_z_PML_2 + kappae_z_PML_2*alphae_z_PML_2) / kappae_z_PML_2;
		}
		ce_z_2(k) = ce_z_2(k) / dz;
		den_ez(k) = 1.0 / (kappae_z_PML_2*dz);
	}
	//  Hz PML2
	double sigh_z_PML_2, alphah_z_PML_2, kappah_z_PML_2;
	for (k = KsMax; k <= Kmax - 1; k++)
	{
		sigh_z_PML_2 = sig_z_max * pow((k + 0.5 - KsMax) / double(nzPML_2), m);
		alphah_z_PML_2 = alpha_z_max*pow((Kmax - k - 0.5) / double(nzPML_2), ma);
		kappah_z_PML_2 = 1.0 + (kappa_z_max - 1.0)*pow((k + 0.5 - KsMax) / double(nzPML_2), m);

		bh_z_2(k) = exp(-(sigh_z_PML_2 / kappah_z_PML_2 + alphah_z_PML_2)*dt / eps0);
		ch_z_2(k) = sigh_z_PML_2*(bh_z_2(k) - 1.0) /
			(sigh_z_PML_2 + kappah_z_PML_2*alphah_z_PML_2) / kappah_z_PML_2;

		ch_z_2(k) = ch_z_2(k) / dz;
		den_hz(k) = 1.0 / (kappah_z_PML_2*dz);
	}
}

void FDTD::update3D_H()
{
	int i, j, k;
	//UPDATE Hx
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax; i++) {       //为什么 i<=Imax
			for (j = Jmin; j <= Jmax - 1; j++) {
				//mm = hxmed(i, j, k);  //计算
				Hx(i, j, k) = DA * Hx(i, j, k) + DB *
					((Ez(i, j, k) - Ez(i, j + 1, k))*den_hy(j) + (Ey(i, j, k + 1) - Ey(i, j, k))*den_hz(k));
			}
		}
	}   //k-loop

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//UPDATE Hy
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax - 1; i++) {
			for (j = Jmin; j <= Jmax; j++) {
				//mm = hymed(i, j, k);
				Hy(i, j, k) = DA * Hy(i, j, k) + DB * ((Ez(i + 1, j, k) - Ez(i, j, k))*den_hx(i) + (Ex(i, j, k) - Ex(i, j, k + 1))*den_hz(k));
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Hz
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax; k++) {
		for (i = Imin; i <= Imax - 1; i++) {
			for (j = Jmin; j <= Jmax - 1; j++) {
				//mm = hzmed(i, j, k);
				Hz(i, j, k) = DA * Hz(i, j, k) + DB * ((Ey(i, j, k) - Ey(i + 1, j, k))*den_hx(i) + (Ex(i, j + 1, k) - Ex(i, j, k))*den_hy(j));
			}
		}
	}

}

void FDTD::update3D_E()
{
	int i, j, k;
	double factor1, factor2;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ex
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin + 1; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax - 1; i++) {
			for (j = Jmin + 1; j <= Jmax - 1; j++) {
				factor1 = CA[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];
				factor2 = CB[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];

				Ex(i, j, k) = factor1 * Ex(i, j, k) + factor2 * ((Hz(i, j, k) - Hz(i, j - 1, k))*den_ey(j) +
					(Hy(i, j, k - 1) - Hy(i, j, k))*den_ez(k));
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ey
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin + 1; k <= Kmax - 1; k++) {
		for (i = Imin + 1; i <= Imax - 1; i++) {
			for (j = Jmin; j <= Jmax - 1; j++) {
				factor1 = CA[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];

				Ey(i, j, k) = factor1 * Ey(i, j, k) + factor2 * ((Hz(i - 1, j, k) - Hz(i, j, k))*den_ex(i) +
					(Hx(i, j, k) - Hx(i, j, k - 1))*den_ez(k));
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ez
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (i = Imin + 1; i <= Imax - 1; i++) {
			for (j = Jmin + 1; j <= Jmax - 1; j++) {
				factor1 = CA[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];

				Ez(i, j, k) = factor1 * Ez(i, j, k) + factor2 * ((Hy(i, j, k) - Hy(i - 1, j, k))*den_ex(i) +
					(Hx(i, j - 1, k) - Hx(i, j, k))*den_ey(j));
			}
		}
	}

}


void FDTD::update3D_CPML_psiH()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Hx
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax; i++) {
			//.....................................................................
			//PML for bottom Hx, j - direction
			//.....................................................................
			for (j = Jmin; j <= JsMin - 1; j++) {
				//mm = hxmed(i, j, k);
				psi_Hxy_1(i, j, k) = bh_y_1(j)*psi_Hxy_1(i, j, k) + ch_y_1(j) *(Ez(i, j, k) - Ez(i, j + 1, k));
				//Hx(i, j, k) = Hx(i, j, k) + DB * psi_Hxy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hx, j - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				//mm = hxmed(i, j, k);
				psi_Hxy_2(i, j, k) = bh_y_2(j)*psi_Hxy_2(i, j, k) + ch_y_2(j) *(Ez(i, j, k) - Ez(i, j + 1, k));
				//Hx(i, j, k) = Hx(i, j, k) + DB * psi_Hxy_2(i, j, k);
			}
		}  // i-loop
	}   //k-loop

	for (i = Imin; i <= Imax; i++) {
		for (j = Jmin; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Hx, k - direction
			//.....................................................................
			for (k = Kmin; k <= KsMin - 1; k++) {
				//mm = hxmed(i, j, k);
				psi_Hxz_1(i, j, k) = bh_z_1(k)*psi_Hxz_1(i, j, k) + ch_z_1(k) *(Ey(i, j, k + 1) - Ey(i, j, k));
				//Hx(i, j, k) = Hx(i, j, k) + DB * psi_Hxz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hx, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				//mm = hxmed(i, j, k);
				psi_Hxz_2(i, j, k) = bh_z_2(k)*psi_Hxz_2(i, j, k) + ch_z_2(k) *(Ey(i, j, k + 1) - Ey(i, j, k));
				//Hx(i, j, k) = Hx(i, j, k) + DB * psi_Hxz_2(i, j, k);
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Hy
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (j = Jmin; j <= Jmax; j++) {
			//.....................................................................
			//PML for bottom Hy, i - direction
			//.....................................................................
			for (i = Imin; i <= IsMin - 1; i++) {
				//mm = hymed(i, j, k);
				psi_Hyx_1(i, j, k) = bh_x_1(i)*psi_Hyx_1(i, j, k) + ch_x_1(i)*(Ez(i + 1, j, k) - Ez(i, j, k));
				//Hy(i, j, k) = Hy(i, j, k) + DB * psi_Hyx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hy, i - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				//mm = hymed(i, j, k);
				psi_Hyx_2(i, j, k) = bh_x_2(i)*psi_Hyx_2(i, j, k) + ch_x_2(i)*(Ez(i + 1, j, k) - Ez(i, j, k));
				//Hy(i, j, k) = Hy(i, j, k) + DB * psi_Hyx_2(i, j, k);
			}
		}
	}
	for (i = Imin; i <= Imax - 1; i++) {
		for (j = Jmin; j <= Jmax; j++) {
			//.....................................................................
			//PML for bottom Hy, k - direction
			//.....................................................................
			for (k = Kmin; k <= KsMin - 1; k++) {
				//mm = hymed(i, j, k);
				psi_Hyz_1(i, j, k) = bh_z_1(k)*psi_Hyz_1(i, j, k) + ch_z_1(k)*(Ex(i, j, k) - Ex(i, j, k + 1));
				//Hy(i, j, k) = Hy(i, j, k) + DB * psi_Hyz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hy, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				//mm = hymed(i, j, k);
				psi_Hyz_2(i, j, k) = bh_z_2(k)*psi_Hyz_2(i, j, k) + ch_z_2(k)*(Ex(i, j, k) - Ex(i, j, k + 1));
				//Hy(i, j, k) = Hy(i, j, k) + DB * psi_Hyz_2(i, j, k);
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Hz
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax; k++) {
		for (j = Jmin; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Hz, x - direction
			//.....................................................................
			for (i = Imin; i <= IsMin - 1; i++) {
				//mm = hzmed(i, j, k);
				psi_Hzx_1(i, j, k) = bh_x_1(i)*psi_Hzx_1(i, j, k) + ch_x_1(i) *(Ey(i, j, k) - Ey(i + 1, j, k));
				//Hz(i, j, k) = Hz(i, j, k) + DB * psi_Hzx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hz, x - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				//mm = hzmed(i, j, k);
				psi_Hzx_2(i, j, k) = bh_x_2(i)*psi_Hzx_2(i, j, k) + ch_x_2(i) *(Ey(i, j, k) - Ey(i + 1, j, k));
				//Hz(i, j, k) = Hz(i, j, k) + DB * psi_Hzx_2(i, j, k);
			}
		}
		for (i = Imin; i <= Imax - 1; i++) {
			//.....................................................................
			//PML for bottom Hz, y - direction
			//.....................................................................
			for (j = Jmin; j <= JsMin - 1; j++) {
				//mm = hzmed(i, j, k);
				psi_Hzy_1(i, j, k) = bh_y_1(j)*psi_Hzy_1(i, j, k) + ch_y_1(j)*(Ex(i, j + 1, k) - Ex(i, j, k));
				//Hz(i, j, k) = Hz(i, j, k) + DB * psi_Hzy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hz, y - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				//mm = hzmed(i, j, k);
				psi_Hzy_2(i, j, k) = bh_y_2(j)*psi_Hzy_2(i, j, k) + ch_y_2(j)*(Ex(i, j + 1, k) - Ex(i, j, k));
				//Hz(i, j, k) = Hz(i, j, k) + DB * psi_Hzy_2(i, j, k);
			}
		}
	}

}

void FDTD::update3D_CPML_psiE()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ex
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin + 1; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax - 1; i++) {
			//.....................................................................
			//PML for bottom Ex, j - direction
			//.....................................................................
			for (j = Jmin + 1; j <= JsMin; j++) {
				//factor2 = CB[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];

				psi_Exy_1(i, j, k) = be_y_1(j)*psi_Exy_1(i, j, k) + ce_y_1(j) *(Hz(i, j, k) - Hz(i, j - 1, k));
				//Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ex, j - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				//factor2 = CB[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];

				psi_Exy_2(i, j, k) = be_y_2(j)*psi_Exy_2(i, j, k) + ce_y_2(j) *(Hz(i, j, k) - Hz(i, (j - 1), k));
				//Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exy_2(i, j, k);
			}
		}
	}
	for (i = Imin; i <= Imax - 1; i++) {
		for (j = Jmin + 1; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Ex, k - direction
			//.....................................................................
			for (k = Kmin + 1; k <= KsMin; k++) {
				//factor2 = CB[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];

				psi_Exz_1(i, j, k) = be_z_1(k)*psi_Exz_1(i, j, k) + ce_z_1(k) *(Hy(i, j, k - 1) - Hy(i, j, k));
				//Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ex, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				//factor2 = CB[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];

				psi_Exz_2(i, j, k) = be_z_2(k)*psi_Exz_2(i, j, k) + ce_z_2(k) *(Hy(i, j, k - 1) - Hy(i, j, k));
				//Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exz_2(i, j, k);
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ey
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin + 1; k <= Kmax - 1; k++) {
		for (j = Jmin; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Ey, i - direction
			//.....................................................................
			for (i = Imin + 1; i <= IsMin; i++) {
				//factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];

				psi_Eyx_1(i, j, k) = be_x_1(i)*psi_Eyx_1(i, j, k) + ce_x_1(i)*(Hz(i - 1, j, k) - Hz(i, j, k));
				//Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ey, i - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				//factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];

				psi_Eyx_2(i, j, k) = be_x_2(i)*psi_Eyx_2(i, j, k) + ce_x_2(i)*(Hz(i - 1, j, k) - Hz(i, j, k));
				//Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyx_2(i, j, k);
			}
		}
	}
	for (i = Imin + 1; i <= Imax - 1; i++) {
		for (j = Jmin; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Ey, k - direction
			//.....................................................................
			for (k = Kmin + 1; k <= KsMin; k++) {
				//factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];

				psi_Eyz_1(i, j, k) = be_z_1(k)*psi_Eyz_1(i, j, k) + ce_z_1(k)*(Hx(i, j, k) - Hx(i, j, k - 1));
				//Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ey, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				//factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];

				psi_Eyz_2(i, j, k) = be_z_2(k)*psi_Eyz_2(i, j, k) + ce_z_2(k)*(Hx(i, j, k) - Hx(i, j, k - 1));
				//Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyz_2(i, j, k);
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ez
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (j = Jmin + 1; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Ez, x - direction
			//.....................................................................
			for (i = Imin + 1; i <= IsMin; i++) {
				//factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];

				psi_Ezx_1(i, j, k) = be_x_1(i)*psi_Ezx_1(i, j, k) + ce_x_1(i) *(Hy(i, j, k) - Hy(i - 1, j, k));
				//Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ez, x - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				//factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];

				psi_Ezx_2(i, j, k) = be_x_2(i)*psi_Ezx_2(i, j, k) + ce_x_2(i) *(Hy(i, j, k) - Hy(i - 1, j, k));
				//Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezx_2(i, j, k);
			}
		}
		for (i = Imin + 1; i <= Imax - 1; i++) {
			//.....................................................................
			//PML for bottom Ez, y - direction
			//.....................................................................
			for (j = Jmin + 1; j <= JsMin; j++) {
				//factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];

				psi_Ezy_1(i, j, k) = be_y_1(j)*psi_Ezy_1(i, j, k) + ce_y_1(j)*(Hx(i, j - 1, k) - Hx(i, j, k));
				//Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ez, y - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				//factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];

				psi_Ezy_2(i, j, k) = be_y_2(j)*psi_Ezy_2(i, j, k) + ce_y_2(j)*(Hx(i, j - 1, k) - Hx(i, j, k));
				//Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezy_2(i, j, k);
			}
		}
	}

}

void FDTD::update3D_CPML_H()
{
	int i, j, k;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Hx
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax; i++) {
			//.....................................................................
			//PML for bottom Hx, j - direction
			//.....................................................................
			for (j = Jmin; j <= JsMin - 1; j++) {
				//mm = hxmed(i, j, k);
				//psi_Hxy_1(i, j, k) = bh_y_1(j)*psi_Hxy_1(i, j, k) + ch_y_1(j) *(Ez(i, j, k) - Ez(i, j + 1, k));
				Hx(i, j, k) = Hx(i, j, k) + DB * psi_Hxy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hx, j - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				//mm = hxmed(i, j, k);
				//psi_Hxy_2(i, j, k) = bh_y_2(j)*psi_Hxy_2(i, j, k) + ch_y_2(j) *(Ez(i, j, k) - Ez(i, j + 1, k));
				Hx(i, j, k) = Hx(i, j, k) + DB * psi_Hxy_2(i, j, k);
			}
		}  // i-loop
	}   //k-loop

	for (i = Imin; i <= Imax; i++) {
		for (j = Jmin; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Hx, k - direction
			//.....................................................................
			for (k = Kmin; k <= KsMin - 1; k++) {
				//mm = hxmed(i, j, k);
				//psi_Hxz_1(i, j, k) = bh_z_1(k)*psi_Hxz_1(i, j, k) + ch_z_1(k) *(Ey(i, j, k + 1) - Ey(i, j, k));
				Hx(i, j, k) = Hx(i, j, k) + DB * psi_Hxz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hx, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				//mm = hxmed(i, j, k);
				//psi_Hxz_2(i, j, k) = bh_z_2(k)*psi_Hxz_2(i, j, k) + ch_z_2(k) *(Ey(i, j, k + 1) - Ey(i, j, k));
				Hx(i, j, k) = Hx(i, j, k) + DB * psi_Hxz_2(i, j, k);
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Hy
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (j = Jmin; j <= Jmax; j++) {
			//.....................................................................
			//PML for bottom Hy, i - direction
			//.....................................................................
			for (i = Imin; i <= IsMin - 1; i++) {
				//mm = hymed(i, j, k);
				//psi_Hyx_1(i, j, k) = bh_x_1(i)*psi_Hyx_1(i, j, k) + ch_x_1(i)*(Ez(i + 1, j, k) - Ez(i, j, k));
				Hy(i, j, k) = Hy(i, j, k) + DB * psi_Hyx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hy, i - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				//mm = hymed(i, j, k);
				//psi_Hyx_2(i, j, k) = bh_x_2(i)*psi_Hyx_2(i, j, k) + ch_x_2(i)*(Ez(i + 1, j, k) - Ez(i, j, k));
				Hy(i, j, k) = Hy(i, j, k) + DB * psi_Hyx_2(i, j, k);
			}
		}
	}
	for (i = Imin; i <= Imax - 1; i++) {
		for (j = Jmin; j <= Jmax; j++) {
			//.....................................................................
			//PML for bottom Hy, k - direction
			//.....................................................................
			for (k = Kmin; k <= KsMin - 1; k++) {
				//mm = hymed(i, j, k);
				//psi_Hyz_1(i, j, k) = bh_z_1(k)*psi_Hyz_1(i, j, k) + ch_z_1(k)*(Ex(i, j, k) - Ex(i, j, k + 1));
				Hy(i, j, k) = Hy(i, j, k) + DB * psi_Hyz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hy, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				//mm = hymed(i, j, k);
				//psi_Hyz_2(i, j, k) = bh_z_2(k)*psi_Hyz_2(i, j, k) + ch_z_2(k)*(Ex(i, j, k) - Ex(i, j, k + 1));
				Hy(i, j, k) = Hy(i, j, k) + DB * psi_Hyz_2(i, j, k);
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Hz
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax; k++) {
		for (j = Jmin; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Hz, x - direction
			//.....................................................................
			for (i = Imin; i <= IsMin - 1; i++) {
				//mm = hzmed(i, j, k);
				//psi_Hzx_1(i, j, k) = bh_x_1(i)*psi_Hzx_1(i, j, k) + ch_x_1(i) *(Ey(i, j, k) - Ey(i + 1, j, k));
				Hz(i, j, k) = Hz(i, j, k) + DB * psi_Hzx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hz, x - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				//mm = hzmed(i, j, k);
				//psi_Hzx_2(i, j, k) = bh_x_2(i)*psi_Hzx_2(i, j, k) + ch_x_2(i) *(Ey(i, j, k) - Ey(i + 1, j, k));
				Hz(i, j, k) = Hz(i, j, k) + DB * psi_Hzx_2(i, j, k);
			}
		}
		for (i = Imin; i <= Imax - 1; i++) {
			//.....................................................................
			//PML for bottom Hz, y - direction
			//.....................................................................
			for (j = Jmin; j <= JsMin - 1; j++) {
				//mm = hzmed(i, j, k);
				//psi_Hzy_1(i, j, k) = bh_y_1(j)*psi_Hzy_1(i, j, k) + ch_y_1(j)*(Ex(i, j + 1, k) - Ex(i, j, k));
				Hz(i, j, k) = Hz(i, j, k) + DB * psi_Hzy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Hz, y - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				//mm = hzmed(i, j, k);
				//psi_Hzy_2(i, j, k) = bh_y_2(j)*psi_Hzy_2(i, j, k) + ch_y_2(j)*(Ex(i, j + 1, k) - Ex(i, j, k));
				Hz(i, j, k) = Hz(i, j, k) + DB * psi_Hzy_2(i, j, k);
			}
		}
	}

}

void FDTD::update3D_CPML_E()
{
	int i, j, k;
	double factor2;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ex
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin + 1; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax - 1; i++) {
			//.....................................................................
			//PML for bottom Ex, j - direction
			//.....................................................................
			for (j = Jmin + 1; j <= JsMin; j++) {
				factor2 = CB[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];

				//psi_Exy_1(i, j, k) = be_y_1(j)*psi_Exy_1(i, j, k) + ce_y_1(j) *(Hz(i, j, k) - Hz(i, j - 1, k));
				Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ex, j - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				factor2 = CB[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];

				//psi_Exy_2(i, j, k) = be_y_2(j)*psi_Exy_2(i, j, k) + ce_y_2(j) *(Hz(i, j, k) - Hz(i, (j - 1), k));
				Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exy_2(i, j, k);
			}
		}
	}
	for (i = Imin; i <= Imax - 1; i++) {
		for (j = Jmin + 1; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Ex, k - direction
			//.....................................................................
			for (k = Kmin + 1; k <= KsMin; k++) {
				factor2 = CB[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];

				//psi_Exz_1(i, j, k) = be_z_1(k)*psi_Exz_1(i, j, k) + ce_z_1(k) *(Hy(i, j, k - 1) - Hy(i, j, k));
				Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ex, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				factor2 = CB[ob(i, j, k)][ob(i, j, k - 1)][ob(i, j - 1, k)][ob(i, j - 1, k - 1)];

				//psi_Exz_2(i, j, k) = be_z_2(k)*psi_Exz_2(i, j, k) + ce_z_2(k) *(Hy(i, j, k - 1) - Hy(i, j, k));
				Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exz_2(i, j, k);
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ey
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin + 1; k <= Kmax - 1; k++) {
		for (j = Jmin; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Ey, i - direction
			//.....................................................................
			for (i = Imin + 1; i <= IsMin; i++) {
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];

				//psi_Eyx_1(i, j, k) = be_x_1(i)*psi_Eyx_1(i, j, k) + ce_x_1(i)*(Hz(i - 1, j, k) - Hz(i, j, k));
				Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ey, i - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];

				//psi_Eyx_2(i, j, k) = be_x_2(i)*psi_Eyx_2(i, j, k) + ce_x_2(i)*(Hz(i - 1, j, k) - Hz(i, j, k));
				Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyx_2(i, j, k);
			}
		}
	}
	for (i = Imin + 1; i <= Imax - 1; i++) {
		for (j = Jmin; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Ey, k - direction
			//.....................................................................
			for (k = Kmin + 1; k <= KsMin; k++) {
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];

				//psi_Eyz_1(i, j, k) = be_z_1(k)*psi_Eyz_1(i, j, k) + ce_z_1(k)*(Hx(i, j, k) - Hx(i, j, k - 1));
				Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ey, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j, k - 1)][ob(i - 1, j, k - 1)];

				//psi_Eyz_2(i, j, k) = be_z_2(k)*psi_Eyz_2(i, j, k) + ce_z_2(k)*(Hx(i, j, k) - Hx(i, j, k - 1));
				Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyz_2(i, j, k);
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ez
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (j = Jmin + 1; j <= Jmax - 1; j++) {
			//.....................................................................
			//PML for bottom Ez, x - direction
			//.....................................................................
			for (i = Imin + 1; i <= IsMin; i++) {
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];

				//psi_Ezx_1(i, j, k) = be_x_1(i)*psi_Ezx_1(i, j, k) + ce_x_1(i) *(Hy(i, j, k) - Hy(i - 1, j, k));
				Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ez, x - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];

				//psi_Ezx_2(i, j, k) = be_x_2(i)*psi_Ezx_2(i, j, k) + ce_x_2(i) *(Hy(i, j, k) - Hy(i - 1, j, k));
				Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezx_2(i, j, k);
			}
		}
		for (i = Imin + 1; i <= Imax - 1; i++) {
			//.....................................................................
			//PML for bottom Ez, y - direction
			//.....................................................................
			for (j = Jmin + 1; j <= JsMin; j++) {
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];

				//psi_Ezy_1(i, j, k) = be_y_1(j)*psi_Ezy_1(i, j, k) + ce_y_1(j)*(Hx(i, j - 1, k) - Hx(i, j, k));
				Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ez, y - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				factor2 = CB[ob(i, j, k)][ob(i - 1, j, k)][ob(i, j - 1, k)][ob(i - 1, j - 1, k)];

				//psi_Ezy_2(i, j, k) = be_y_2(j)*psi_Ezy_2(i, j, k) + ce_y_2(j)*(Hx(i, j - 1, k) - Hx(i, j, k));
				Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezy_2(i, j, k);
			}
		}
	}

}