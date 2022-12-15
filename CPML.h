#ifndef __CPML_H__
#define __CPML_H__
//////////////////////////////////////////////////////////////////////
//PURPOSE:
//	CPML absorbing boundary conditions for 3 - D FDTD code.
//NOTES:
//Reference:
//	[1]  葛德彪, 闫玉波. (2011). 电磁波时域有限差分方法（第三版）.西安电子科技大学出版社.
//  [2]  Taflove, A., &Hagness, S.C. (2005).Computational electrodynamics 
//       : the finite - difference time - domain method 3rd ed. Artech house.
//Record of revisions:
//   1. 06 May 2020 by Yang chao, Email: yangchaophy@foxmail.com
//   2. 02 November 2022 by Yang chao, Email: yangchaophy@foxmail.com
///////////////////////////////////////////////////////////////////////////

#include<vector>
#include "matrix.h"
#include"physical_constants.h"


struct CPMLPrivate
{
	//Specify the CPML Order and Other Parameters
	double sig_factor = 1.0; //取值范围 sig_factor=0.7 ~ 1.5
	int m = 3;               //取值范围 m=2,3,4
	int ma = 1;
	double sig_max = sig_factor*0.8*(m + 1) / Z0;
	//[1]P107  alpha=0~0.05 
	double alpha_max = 0.0;
	//[1]P107  kappa=5~11
	double kappa_max = 1.0;
};


class CPML
{
private:
	int nxPML_1, nxPML_2, nyPML_1, nyPML_2, nzPML_1, nzPML_2;
	int IsMin, IsMax, JsMin, JsMax, KsMin, KsMax;  //散射场区域边界
	int Imin, Imax, Jmin, Jmax, Kmin, Kmax;  //  总区域
	
	CPMLPrivate CPMLpar;

public:
	Matrix<double> psi_Exy_1, psi_Exy_2;
	Matrix<double> psi_Exz_1, psi_Exz_2;
	Matrix<double> psi_Eyx_1, psi_Eyx_2;
	Matrix<double> psi_Eyz_1, psi_Eyz_2;
	Matrix<double> psi_Ezx_1, psi_Ezx_2;
	Matrix<double> psi_Ezy_1, psi_Ezy_2;

	Matrix<double> psi_Hxy_1, psi_Hxy_2;
	Matrix<double> psi_Hxz_1, psi_Hxz_2;
	Matrix<double> psi_Hyx_1, psi_Hyx_2;
	Matrix<double> psi_Hyz_1, psi_Hyz_2;
	Matrix<double> psi_Hzx_1, psi_Hzx_2;
	Matrix<double> psi_Hzy_1, psi_Hzy_2;

	Matrix<double> be_x_1, be_x_2;
	Matrix<double> be_y_1, be_y_2;
	Matrix<double> be_z_1, be_z_2;
	Matrix<double> ce_x_1, ce_x_2;
	Matrix<double> ce_y_1, ce_y_2;
	Matrix<double> ce_z_1, ce_z_2;

	Matrix<double> bh_x_1, bh_x_2;
	Matrix<double> bh_y_1, bh_y_2;
	Matrix<double> bh_z_1, bh_z_2;
	Matrix<double> ch_x_1, ch_x_2;
	Matrix<double> ch_y_1, ch_y_2;
	Matrix<double> ch_z_1, ch_z_2;
	
public:
	CPML();
	CPML(int cpmlWidth, int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_);
	CPML(int width_xn, int width_xp, int width_yn, int width_yp, int width_zn, int width_zp, 
		int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_);
	
	void initializeCPML();

	void initCoefficientsCPML(const std::vector<double> epsr, double dx, double dy, double dz, double dt,
		Matrix<double> &den_ex, Matrix<double> &den_ey, Matrix<double> &den_ez,
		Matrix<double> &den_hx, Matrix<double> &den_hy, Matrix<double> &den_hz);
	
	void update3D_CPML_psiE(const Matrix<double> &Hx, const Matrix<double> &Hy, const Matrix<double> &Hz);

	void update3D_CPML_psiH(const Matrix<double> &Ex, const Matrix<double> &Ey, const Matrix<double> &Ez);

	void update3D_CPML_E(const Matrix<double> &CB, const Matrix<int> &ob, Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez);

	void update3D_CPML_H(double DB, Matrix<double> &Hx, Matrix<double>& Hy, Matrix<double> &Hz);


	
private:
	void setCPMLRegion(int width_xn, int width_xp, int width_yn, int width_yp, int width_zn, int width_zp,
		int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_);
};

#endif 