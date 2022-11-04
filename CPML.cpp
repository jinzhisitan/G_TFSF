#include<iostream>
#include<cmath>
#include "CPML.h"

CPML::CPML() :nxPML_1(0), nxPML_2(0), nyPML_1(0), nyPML_2(0), nzPML_1(0), nzPML_2(0)
{
	std::cout << "class CMPL: void constructor\n";
}

CPML::CPML(int cpmlWidth, int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_)
{
	setCPMLRegion(cpmlWidth, cpmlWidth, cpmlWidth, cpmlWidth, cpmlWidth, cpmlWidth,
		IsMin_, IsMax_, JsMin_, JsMax_, KsMin_, KsMax_);
}

CPML::CPML(int width_xn, int width_xp, int width_yn, int width_yp, int width_zn, int width_zp,
	int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_)
{
	setCPMLRegion(width_xn, width_xp, width_yn, width_yp, width_zn, width_zp, IsMin_, IsMax_, JsMin_, JsMax_, KsMin_, KsMax_);
}


void CPML::setCPMLRegion(int width_xn, int width_xp, int width_yn, int width_yp, int width_zn, int width_zp,
	int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_)
{
	nxPML_1 = width_xn;
	nxPML_2 = width_xp;
	nyPML_1 = width_yn;
	nyPML_2 = width_yp;
	nzPML_1 = width_zn;
	nzPML_2 = width_zp;

	IsMin = IsMin_;
	IsMax = IsMax_;
	JsMin = JsMin_;
	JsMax = JsMax_;
	KsMin = KsMin_;
	KsMax = KsMax_;

	Imin = IsMin_ - width_xn;
	Imax = IsMax_ + width_xp;
	Jmin = JsMin_ - width_yn;
	Jmax = JsMax_ + width_yp;
	Kmin = KsMin_ - width_zn;
	Kmax = KsMax_ + width_zp;
}

void CPML::initializeCPML()
{
	psi_Exy_1.Allocate(Imin, Imax - 1, Jmin + 1, JsMin, Kmin + 1, Kmax - 1); psi_Exy_2.Allocate(Imin, Imax - 1, JsMax, Jmax - 1, Kmin + 1, Kmax - 1);
	psi_Exz_1.Allocate(Imin, Imax - 1, Jmin + 1, Jmax - 1, Kmin + 1, KsMin), psi_Exz_2.Allocate(Imin, Imax - 1, Jmin + 1, Jmax - 1, KsMax, Kmax - 1);
	psi_Eyx_1.Allocate(Imin + 1, IsMin, Jmin, Jmax - 1, Kmin + 1, Kmax - 1), psi_Eyx_2.Allocate(IsMax, Imax - 1, Jmin, Jmax - 1, Kmin + 1, Kmax - 1);
	psi_Eyz_1.Allocate(Imin + 1, Imax - 1, Jmin, Jmax - 1, Kmin + 1, KsMin), psi_Eyz_2.Allocate(Imin + 1, Imax - 1, Jmin, Jmax - 1, KsMax, Kmax - 1);
	psi_Ezx_1.Allocate(Imin + 1, IsMin, Jmin + 1, Jmax - 1, Kmin, Kmax - 1), psi_Ezx_2.Allocate(IsMax, Imax - 1, Jmin + 1, Jmax - 1, Kmin, Kmax - 1);
	psi_Ezy_1.Allocate(Imin + 1, Imax - 1, Jmin + 1, JsMin, Kmin, Kmax - 1), psi_Ezy_2.Allocate(Imin + 1, Imax - 1, JsMax, Jmax - 1, Kmin, Kmax - 1);

	psi_Hxy_1.Allocate(Imin, Imax, Jmin, JsMin - 1, Kmin, Kmax - 1), psi_Hxy_2.Allocate(Imin, Imax, JsMax, Jmax - 1, Kmin, Kmax - 1);
	psi_Hxz_1.Allocate(Imin, Imax, Jmin, Jmax - 1, Kmin, KsMin - 1), psi_Hxz_2.Allocate(Imin, Imax, Jmin, Jmax - 1, KsMax, Kmax - 1);
	psi_Hyx_1.Allocate(Imin, IsMin - 1, Jmin, Jmax, Kmin, Kmax - 1), psi_Hyx_2.Allocate(IsMax, Imax - 1, Jmin, Jmax, Kmin, Kmax - 1);
	psi_Hyz_1.Allocate(Imin, Imax - 1, Jmin, Jmax, Kmin, KsMin - 1), psi_Hyz_2.Allocate(Imin, Imax - 1, Jmin, Jmax, KsMax, Kmax - 1);
	psi_Hzx_1.Allocate(Imin, IsMin - 1, Jmin, Jmax - 1, Kmin, Kmax), psi_Hzx_2.Allocate(IsMax, Imax - 1, Jmin, Jmax - 1, Kmin, Kmax);
	psi_Hzy_1.Allocate(Imin, Imax - 1, Jmin, JsMin - 1, Kmin, Kmax), psi_Hzy_2.Allocate(Imin, Imax - 1, JsMax, Jmax - 1, Kmin, Kmax);

	be_x_1.Allocate(Imin + 1, IsMin), ce_x_1.Allocate(Imin + 1, IsMin);
	be_x_2.Allocate(IsMax, Imax - 1), ce_x_2.Allocate(IsMax, Imax - 1);
	be_y_1.Allocate(Jmin + 1, JsMin), ce_y_1.Allocate(Jmin + 1, JsMin);
	be_y_2.Allocate(JsMax, Jmax - 1), ce_y_2.Allocate(JsMax, Jmax - 1);
	be_z_1.Allocate(Kmin + 1, KsMin), ce_z_1.Allocate(Kmin + 1, KsMin);
	be_z_2.Allocate(KsMax, Kmax - 1), ce_z_2.Allocate(KsMax, Kmax - 1);

	bh_x_1.Allocate(Imin, IsMin - 1), ch_x_1.Allocate(Imin, IsMin - 1);
	bh_x_2.Allocate(IsMax, Imax - 1), ch_x_2.Allocate(IsMax, Imax - 1);
	bh_y_1.Allocate(Jmin, JsMin - 1), ch_y_1.Allocate(Jmin, JsMin - 1);
	bh_y_2.Allocate(JsMax, Jmax - 1), ch_y_2.Allocate(JsMax, Jmax - 1);
	bh_z_1.Allocate(Kmin, KsMin - 1), ch_z_1.Allocate(Kmin, KsMin - 1);
	bh_z_2.Allocate(KsMax, Kmax - 1), ch_z_2.Allocate(KsMax, Kmax - 1);
}


void CPML::initCoefficientsCPML(const std::vector<double> epsr, double dx, double dy, double dz, double dt,
	Matrix<double> &den_ex, Matrix<double> &den_ey, Matrix<double> &den_ez,
	Matrix<double> &den_hx, Matrix<double> &den_hy, Matrix<double> &den_hz)
{
	/////////////////////////////////////////////////////////////////////////////
	//Specify the CPML Order and Other Parameters
	//[1]Eq(5-4-50), [2]Eq(7.66),Eq(7.67)
	double sig_x_max = CPMLpar.sig_max / (dx*sqrt(epsr[0]));
	double sig_y_max = CPMLpar.sig_max / (dy*sqrt(epsr[0]));
	double sig_z_max = CPMLpar.sig_max / (dz*sqrt(epsr[0]));
	//[1]P107  alpha=0~0.05 
	double alpha_x_max = CPMLpar.alpha_max, alpha_y_max = alpha_x_max, alpha_z_max = alpha_x_max;
	//[1]P107  kappa=5~11
	double kappa_x_max = CPMLpar.kappa_max, kappa_y_max = kappa_x_max, kappa_z_max = kappa_x_max;
	///////////////////////////////////////////////////////////////////////////////

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//SET CPML PARAMETERS IN EACH DIRECTION
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int i, j, k;
	//  Ex  PML1
	double sige_x_PML_1, alphae_x_PML_1, kappae_x_PML_1;
	for (i = Imin+1; i <= IsMin; i++) {
		//[1]Eq(5-4-49)   [2]Eq(7.60a)
		sige_x_PML_1 = sig_x_max * pow((IsMin - i) / double(nxPML_1), CPMLpar.m);
		//[1]Eq(5-4-52)   [2]Eq(7.79)  
		alphae_x_PML_1 = alpha_x_max*pow((i - Imin) / double(nxPML_1), CPMLpar.ma);     //ÐÞ¸Ä
																				//[1]Eq(5-4-51)   [2]Eq(7.60b)
		kappae_x_PML_1 = 1.0 + (kappa_x_max - 1.0)*pow((IsMin - i) / double(nxPML_1), CPMLpar.m);

		be_x_1(i) = exp(-(sige_x_PML_1 / kappae_x_PML_1 + alphae_x_PML_1)  *dt / eps0);
		if ((sige_x_PML_1 == 0.0) && (alphae_x_PML_1 == 0.0)&& (i == IsMin)) {
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
		sigh_x_PML_1 = sig_x_max * pow((IsMin - i - 0.5) / double(nxPML_1), CPMLpar.m);
		alphah_x_PML_1 = alpha_x_max*pow((i + 0.5 - Imin) / double(nxPML_1), CPMLpar.ma);
		kappah_x_PML_1 = 1.0 + (kappa_x_max - 1.0)* pow((IsMin - i - 0.5) / double(nxPML_1), CPMLpar.m);

		bh_x_1(i) = exp(-(sigh_x_PML_1 / kappah_x_PML_1 + alphah_x_PML_1)*dt / eps0);
		ch_x_1(i) = sigh_x_PML_1*(bh_x_1(i) - 1.0) /
			(sigh_x_PML_1 + kappah_x_PML_1*alphah_x_PML_1) / kappah_x_PML_1;

		ch_x_1(i) = ch_x_1(i) / dx;
		den_hx(i) = 1.0 / (kappah_x_PML_1*dx);
	}


	//  Ex  PML2
	double sige_x_PML_2, alphae_x_PML_2, kappae_x_PML_2;
	for (i = IsMax; i <= Imax-1; i++)
	{
		sige_x_PML_2 = sig_x_max * pow((i - IsMax) / double(nxPML_2), CPMLpar.m);
		alphae_x_PML_2 = alpha_x_max*pow((Imax - i) / double(nxPML_2), CPMLpar.ma);
		kappae_x_PML_2 = 1.0 + (kappa_x_max - 1.0)*pow((i - IsMax) / double(nxPML_2), CPMLpar.m);

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
		sigh_x_PML_2 = sig_x_max * pow((i + 0.5 - IsMax) / double(nxPML_2), CPMLpar.m);
		alphah_x_PML_2 = alpha_x_max*pow((Imax - i - 0.5) / double(nxPML_2), CPMLpar.ma);
		kappah_x_PML_2 = 1.0 + (kappa_x_max - 1.0)*pow((i + 0.5 - IsMax) / double(nxPML_2), CPMLpar.m);

		bh_x_2(i) = exp(-(sigh_x_PML_2 / kappah_x_PML_2 + alphah_x_PML_2)*dt / eps0);
		ch_x_2(i) = sigh_x_PML_2*(bh_x_2(i) - 1.0) /
			(sigh_x_PML_2 + kappah_x_PML_2*alphah_x_PML_2) / kappah_x_PML_2;

		ch_x_2(i) = ch_x_2(i) / dx;
		den_hx(i) = 1.0 / (kappah_x_PML_2*dx);
	}

	//  Ey  PML1  
	double sige_y_PML_1, alphae_y_PML_1, kappae_y_PML_1;
	for (j = Jmin+1; j <= JsMin; j++)
	{
		sige_y_PML_1 = sig_y_max * pow((JsMin - j) / double(nyPML_1), CPMLpar.m);
		alphae_y_PML_1 = alpha_y_max*pow((j - Jmin) / double(nyPML_1), CPMLpar.ma);
		kappae_y_PML_1 = 1.0 + (kappa_y_max - 1.0)*pow((JsMin - j) / double(nyPML_1), CPMLpar.m);

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
		sigh_y_PML_1 = sig_y_max * pow((JsMin - j - 0.5) / double(nyPML_1), CPMLpar.m);
		alphah_y_PML_1 = alpha_y_max*pow((j + 0.5 - Jmin) / double(nyPML_1), CPMLpar.ma);
		kappah_y_PML_1 = 1.0 + (kappa_y_max - 1.0)*pow((JsMin - j - 0.5) / double(nyPML_1), CPMLpar.m);

		bh_y_1(j) = exp(-(sigh_y_PML_1 / kappah_y_PML_1 + alphah_y_PML_1)*dt / eps0);
		ch_y_1(j) = sigh_y_PML_1*(bh_y_1(j) - 1.0) /
			(sigh_y_PML_1 + kappah_y_PML_1*alphah_y_PML_1) / kappah_y_PML_1;

		ch_y_1(j) = ch_y_1(j) / dy;
		den_hy(j) = 1.0 / (kappah_y_PML_1*dy);
	}

	//  Ey PML2
	double sige_y_PML_2, alphae_y_PML_2, kappae_y_PML_2;
	for (j = JsMax; j <= Jmax-1; j++)
	{
		sige_y_PML_2 = sig_y_max * pow((j - JsMax) / double(nyPML_2), CPMLpar.m);
		alphae_y_PML_2 = alpha_y_max*pow((Jmax - j) / double(nyPML_2), CPMLpar.ma);
		kappae_y_PML_2 = 1.0 + (kappa_y_max - 1.0)*pow((j - JsMax) / double(nyPML_2), CPMLpar.m);

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
		sigh_y_PML_2 = sig_y_max * pow((j + 0.5 - JsMax) / double(nyPML_2), CPMLpar.m);
		alphah_y_PML_2 = alpha_y_max*pow((Jmax - j - 0.5) / double(nyPML_2), CPMLpar.ma);
		kappah_y_PML_2 = 1.0 + (kappa_y_max - 1.0)*pow((j + 0.5 - JsMax) / double(nyPML_2), CPMLpar.m);

		bh_y_2(j) = exp(-(sigh_y_PML_2 / kappah_y_PML_2 + alphah_y_PML_2)*dt / eps0);
		ch_y_2(j) = sigh_y_PML_2*(bh_y_2(j) - 1.0) /
			(sigh_y_PML_2 + kappah_y_PML_2*alphah_y_PML_2) / kappah_y_PML_2;

		ch_y_2(j) = ch_y_2(j) / dy;
		den_hy(j) = 1.0 / (kappah_y_PML_2*dy);
	}

	//  Ez PML1
	double sige_z_PML_1, alphae_z_PML_1, kappae_z_PML_1;
	for (k = Kmin+1; k <= KsMin; k++) {
		sige_z_PML_1 = sig_z_max * pow((KsMin - k) / double(nzPML_1), CPMLpar.m);
		alphae_z_PML_1 = alpha_z_max*pow((k - Kmin) / double(nzPML_1), CPMLpar.ma);
		kappae_z_PML_1 = 1.0 + (kappa_z_max - 1.0)*pow((KsMin - k) / double(nzPML_1), CPMLpar.m);

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
		sigh_z_PML_1 = sig_z_max * pow((KsMin - k - 0.5) / double(nzPML_1), CPMLpar.m);
		alphah_z_PML_1 = alpha_z_max*pow((k + 0.5 - Kmin) / double(nzPML_1), CPMLpar.ma);
		kappah_z_PML_1 = 1.0 + (kappa_z_max - 1.0)*pow((KsMin - k - 0.5) / double(nzPML_1), CPMLpar.m);

		bh_z_1(k) = exp(-(sigh_z_PML_1 / kappah_z_PML_1 + alphah_z_PML_1)*dt / eps0);
		ch_z_1(k) = sigh_z_PML_1*(bh_z_1(k) - 1.0) /
			(sigh_z_PML_1 + kappah_z_PML_1*alphah_z_PML_1) / kappah_z_PML_1;

		ch_z_1(k) = ch_z_1(k) / dz;
		den_hz(k) = 1.0 / (kappah_z_PML_1*dz);
	}

	//  Ez PML2
	double sige_z_PML_2, alphae_z_PML_2, kappae_z_PML_2;
	for (k = KsMax; k <= Kmax-1; k++)
	{
		sige_z_PML_2 = sig_z_max * pow((k - KsMax) / double(nzPML_2), CPMLpar.m);
		alphae_z_PML_2 = alpha_z_max*pow((Kmax - k) / double(nzPML_2), CPMLpar.ma);
		kappae_z_PML_2 = 1.0 + (kappa_z_max - 1.0)*pow((k - KsMax) / double(nzPML_2), CPMLpar.m);

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
		sigh_z_PML_2 = sig_z_max * pow((k + 0.5 - KsMax) / double(nzPML_2), CPMLpar.m);
		alphah_z_PML_2 = alpha_z_max*pow((Kmax - k - 0.5) / double(nzPML_2), CPMLpar.ma);
		kappah_z_PML_2 = 1.0 + (kappa_z_max - 1.0)*pow((k + 0.5 - KsMax) / double(nzPML_2), CPMLpar.m);

		bh_z_2(k) = exp(-(sigh_z_PML_2 / kappah_z_PML_2 + alphah_z_PML_2)*dt / eps0);
		ch_z_2(k) = sigh_z_PML_2*(bh_z_2(k) - 1.0) /
			(sigh_z_PML_2 + kappah_z_PML_2*alphah_z_PML_2) / kappah_z_PML_2;

		ch_z_2(k) = ch_z_2(k) / dz;
		den_hz(k) = 1.0 / (kappah_z_PML_2*dz);
	}
}


void CPML::update3D_CPML_psiE(const Matrix<double> &Hx, const Matrix<double> &Hy, const Matrix<double> &Hz)
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

void CPML::update3D_CPML_psiH(const Matrix<double> &Ex, const Matrix<double> &Ey, const Matrix<double> &Ez)
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


void CPML::update3D_CPML_E(const Matrix<double> &CB, const Matrix<int> &ob, 
	Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez)
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
				factor2 = CB(ob(i, j, k), ob(i, j, k - 1), ob(i, j - 1, k), ob(i, j - 1, k - 1));

				//psi_Exy_1(i, j, k) = be_y_1(j)*psi_Exy_1(i, j, k) + ce_y_1(j) *(Hz(i, j, k) - Hz(i, j - 1, k));
				Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ex, j - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				factor2 = CB(ob(i, j, k), ob(i, j, k - 1), ob(i, j - 1, k), ob(i, j - 1, k - 1));

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
				factor2 = CB(ob(i, j, k), ob(i, j, k - 1), ob(i, j - 1, k), ob(i, j - 1, k - 1));

				//psi_Exz_1(i, j, k) = be_z_1(k)*psi_Exz_1(i, j, k) + ce_z_1(k) *(Hy(i, j, k - 1) - Hy(i, j, k));
				Ex(i, j, k) = Ex(i, j, k) + factor2 * psi_Exz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ex, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				factor2 = CB(ob(i, j, k), ob(i, j, k - 1), ob(i, j - 1, k), ob(i, j - 1, k - 1));

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
				factor2 = CB(ob(i, j, k), ob(i - 1, j, k), ob(i, j, k - 1), ob(i - 1, j, k - 1));

				//psi_Eyx_1(i, j, k) = be_x_1(i)*psi_Eyx_1(i, j, k) + ce_x_1(i)*(Hz(i - 1, j, k) - Hz(i, j, k));
				Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ey, i - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				factor2 = CB(ob(i, j, k), ob(i - 1, j, k), ob(i, j, k - 1), ob(i - 1, j, k - 1));

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
				factor2 = CB(ob(i, j, k), ob(i - 1, j, k), ob(i, j, k - 1), ob(i - 1, j, k - 1));

				//psi_Eyz_1(i, j, k) = be_z_1(k)*psi_Eyz_1(i, j, k) + ce_z_1(k)*(Hx(i, j, k) - Hx(i, j, k - 1));
				Ey(i, j, k) = Ey(i, j, k) + factor2 * psi_Eyz_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ey, k - direction
			//.....................................................................
			for (k = KsMax; k <= Kmax - 1; k++) {
				factor2 = CB(ob(i, j, k), ob(i - 1, j, k), ob(i, j, k - 1), ob(i - 1, j, k - 1));

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
				factor2 = CB(ob(i, j, k), ob(i - 1, j, k), ob(i, j - 1, k), ob(i - 1, j - 1, k));

				//psi_Ezx_1(i, j, k) = be_x_1(i)*psi_Ezx_1(i, j, k) + ce_x_1(i) *(Hy(i, j, k) - Hy(i - 1, j, k));
				Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezx_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ez, x - direction
			//.....................................................................
			for (i = IsMax; i <= Imax - 1; i++) {
				factor2 = CB(ob(i, j, k), ob(i - 1, j, k), ob(i, j - 1, k), ob(i - 1, j - 1, k));

				//psi_Ezx_2(i, j, k) = be_x_2(i)*psi_Ezx_2(i, j, k) + ce_x_2(i) *(Hy(i, j, k) - Hy(i - 1, j, k));
				Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezx_2(i, j, k);
			}
		}
		for (i = Imin + 1; i <= Imax - 1; i++) {
			//.....................................................................
			//PML for bottom Ez, y - direction
			//.....................................................................
			for (j = Jmin + 1; j <= JsMin; j++) {
				factor2 = CB(ob(i, j, k), ob(i - 1, j, k), ob(i, j - 1, k), ob(i - 1, j - 1, k));

				//psi_Ezy_1(i, j, k) = be_y_1(j)*psi_Ezy_1(i, j, k) + ce_y_1(j)*(Hx(i, j - 1, k) - Hx(i, j, k));
				Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezy_1(i, j, k);
			}
			//.....................................................................
			//PML for top Ez, y - direction
			//.....................................................................
			for (j = JsMax; j <= Jmax - 1; j++) {
				factor2 = CB(ob(i, j, k), ob(i - 1, j, k), ob(i, j - 1, k), ob(i - 1, j - 1, k));

				//psi_Ezy_2(i, j, k) = be_y_2(j)*psi_Ezy_2(i, j, k) + ce_y_2(j)*(Hx(i, j - 1, k) - Hx(i, j, k));
				Ez(i, j, k) = Ez(i, j, k) + factor2 * psi_Ezy_2(i, j, k);
			}
		}
	}
}

void CPML::update3D_CPML_H(double DB, Matrix<double> &Hx, Matrix<double>& Hy, Matrix<double> &Hz)
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