#ifndef __FDTD_H__
#define __FDTD_H__
/////////////////////////////////////////////////////////////////////////////////////////////
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
//   2. 02 November 2022 by Yang chao, Email: yangchaophy@foxmail.com
//////////////////////////////////////////////////////////////////////////////////////////////

#include<vector>
#include "matrix.h"
#include "CPML.h"
#include "TFSF.h"

class FDTD
{
private:
	//Specify Number of Time Steps and Grid Size Parameters
	int IsMin, IsMax, JsMin, JsMax, KsMin, KsMax;
	int Imin, Imax, Jmin, Jmax, Kmin, Kmax;  //  总区域
	int tmax; //total number of time steps

	double dx, dy, dz, dt; 
	
	std::vector<double> epsr;
	std::vector<double> sig;

	Matrix<double> Ex, Ey, Ez;
	Matrix<double> Hx, Hy, Hz;
	//denominators for the update equations
	Matrix<double> den_ex, den_ey, den_ez;  
	Matrix<double> den_hx, den_hy, den_hz;
	
	//网格点(i+0.5,j+0.5,k+0.5)点介质的编号。
	Matrix<int> ob;  

	Matrix<double> CA, CB;
	double DA, DB;
	
	TFSF Einc;
	CPML pml;
	
public:
	FDTD();
	FDTD(int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_,
		int xPML1, int xPML2, int yPML1, int yPML2, int zPML1, int zPML2, int tmax_,
		double dx_, double dy_, double dz_, double dt_, std::vector<double> epsr_, std::vector<double> sig_,
		double alpha, double thi, double phi, int ItMin, int ItMax, int JtMin, int JtMax, int KtMin, int KtMax,int IncStart, int IncEnd);
	
	void setUp();
	void compute(); //E & H Field update equation

private:
	void initialize();  //Memory initialization
	void initCoefficients(); 
	void update3D_E();
	void update3D_H();

	void IsMedia();
};
#endif 