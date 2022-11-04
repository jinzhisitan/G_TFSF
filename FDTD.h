#pragma once
//PURPOSE:
//	3 - D FDTD code with CPML absorbing boundary conditions.  
//	This class implements the finite-difference time-domain solution of Maxwell's curl equations 
//  over a three - dimensional Cartesian space lattice.The grid is terminated by CPML absorbing
//  boundary conditions.
//NOTES:
//	1. 需要调用class CPML计算PML内的电磁场，调用class TSFS加载入射平面波。
//Reference:
//	[1]  葛德彪, 闫玉波. (2011). 电磁波时域有限差分方法（第三版）.西安电子科技大学出版社.
//  [2]  Taflove, A., &Hagness, S.C. (2005).Computational electrodynamics 
//       : the finite - difference time - domain method 3rd ed. Artech house.
//Record of revisions:
//   1. 06 May 2020 by Yang chao, Email: yangchaophy@foxmail.com
//      Northwest Institute of Nuclear Technology, Xi'an; 710024, China

#include <iostream>
#include "TFSF.h"
#include "matrix.h"
#include "physical_constants.h"

using namespace std;
class FDTD
{
private:
	//Specify Number of Time Steps and Grid Size Parameters
	//散射场区域边界
	int IsMin, IsMax, JsMin, JsMax, KsMin, KsMax;  
	//PML thickness in each direction
	int nxPML_1, nxPML_2, nyPML_1, nyPML_2, nzPML_1, nzPML_2;  
	int Imin, Jmin, Kmin, Imax, Jmax, Kmax;  //  总区域
	int tmax; //total number of time steps
	double dx, dy, dz, dt; 
	double epsr[MediaNo];
	double sig[MediaNo];

	Matrix<double> Ex, Ey, Ez;
	Matrix<double> Hx, Hy, Hz;
	//denominators for the update equations
	Matrix<double> den_ex, den_ey, den_ez;  
	Matrix<double> den_hx, den_hy, den_hz;
	
	//网格点(i+0.5,j+0.5,k+0.5)点介质的编号。
	Matrix<int> ob;  

	double CA[MediaNo][MediaNo][MediaNo][MediaNo], CB[MediaNo][MediaNo][MediaNo][MediaNo];
	double DA, DB;

/////////////////////////////////////////////////////////////
/////      CPML
//////////////////////////////////////////////////////////////
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

	TFSF Einc;

public:
	FDTD();
	FDTD(int tmax_, int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_, int xPML1, int xPML2,
		int yPML1,int yPML2, int zPML1,int zPML2, double dt_, double dx_, double dy_, double dz_, double epsr_[], double sig_[],
		double alpha, double thi, double phi, int IncStart, int IncEnd, int ItMin, int ItMax, int JtMin, int JtMax, int KtMin, int KtMax);
	
////////////////////////////////////////////////////
////  初始化
	void initCoeficients();
	void IsMedia();

/////////////////////////////////////////////////////
////  计算
	void compute(); //E & H Field update equation
	
	void update3D_E();
	void update3D_H();
	void update3D_CPML_psiE();
	void update3D_CPML_psiH();
	void update3D_CPML_E();
	void update3D_CPML_H();
};
