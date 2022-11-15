//////////////////////////////////////////////////////////////////////////////////////////////////
//PURPOSE:
//	 3 - D FDTD code with CPML absorbing boundary conditions
//   This code implements the finite-difference time-domain solution of Maxwell's curl equations 
//   over a three-dimensional Cartesian space lattice. The grid is terminated by CPML absorbing
//   boundary conditions. The number of grid cells as well as the thickness of
//   the PML in each Cartesian direction can be varied independently.
//NOTES:
//	1. 网格的编号在FDTD_compute.cpp中的函数IsMedia设置。
//	2. 入射波在TFSF_source.cpp中的函数source设置。
//Reference:
//	[1]  葛德彪, 闫玉波. (2011). 电磁波时域有限差分方法（第三版）.西安电子科技大学出版社.
//  [2]  Taflove, A., &Hagness, S.C. (2005).Computational electrodynamics 
//       : the finite - difference time - domain method 3rd ed. Artech house.
//Record of revisions:
//   1. 06 June 2020 by Yang chao, Email: yangchaophy@foxmail.com
//      Northwest Institute of Nuclear Technology, Xi'an; 710024, China
//	 2. 19 August 2020 by Yang chao
//      采用相邻网格电磁参数的等效值（平均值）。ref [1]。P21:section 2.5, P175,section 8.4.1。
//	 3. 12 April 2022 by Yang chao
//      修正一维辅助平面波加载的网格间隔，使之与三维FDTD色散相同
//////////////////////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<ctime>
#include"physical_constants.h"
#include "FDTD.h"

using namespace std;

int main()
{
	clock_t startTime = clock();

	//Specify Grid Cell Size in Each Direction and Calculate the
	//Resulting Courant - Stable Time Step
	//cell size in each direction
	double dx, dy, dz, dt;
	dx = 0.025, dy = dx, dz = dx;
	//time step increment
	//dt = 0.99 / (c0*sqrt(1.0 / pow(dx, 2) + 1.0 / pow(dy, 2) + 1.0 / pow(dz, 2)));    //Eq(3-2-8)	
	dt = dx / 2.0 / c0;
	
	int IsMin, IsMax, JsMin, JsMax, KsMin, KsMax;  //散射场区域边界
	IsMin = -40, IsMax = 40;
	JsMin = IsMin, JsMax = IsMax;
	KsMin =-10, KsMax = 45;
	int tmax=300; //total number of time steps

	//Specify the CPML Thickness in Each Direction(Value of Zero
	//Corresponds to No PML, and the Grid is Terminated with a PEC)
	//PML thickness in each direction
	int nxPML_1, nxPML_2, nyPML_1, nyPML_2, nzPML_1, nzPML_2;
	nxPML_1 = 16, nxPML_2 = nxPML_1, nyPML_1 = nxPML_1, nyPML_2 = nxPML_1, nzPML_1 = nxPML_1, nzPML_2 = nxPML_1;
	
	//Specify Material Relative Permittivity and Conductivity
	vector<double> epsr { 1.0, 10.0 };
	vector<double> sig { 0.0, 0.001 };

	//************************************************************************************
	//平面波加载
	//************************************************************************************
	int TFSFGrid = -12;
	int ItMin = IsMin +TFSFGrid;   //connective boundary
	int ItMax = IsMax -TFSFGrid;
	int JtMin = JsMin +TFSFGrid;
	int JtMax = JsMax -TFSFGrid;
	int KtMin = KsMin -14;
	int KtMax = KsMax -5;
	
	double alpha = 0.0 / 180.0*pi;
	double thi = 180.0 / 180.0*pi;
	double phi = 180.0 / 180.0*pi;
	int IncStart = -500;
	int IncEnd = 500;
	//************************************************************************************

	FDTD fdtd(IsMin, IsMax, JsMin, JsMax, KsMin, KsMax, nxPML_1, nxPML_2, nyPML_1, nyPML_2, nzPML_1, nzPML_2, tmax,
		dx, dy, dz, dt, epsr, sig, alpha, thi, phi, ItMin, ItMax, JtMin, JtMax, KtMin, KtMax, IncStart, IncEnd);
	

	fdtd.setUp();
	fdtd.compute();
	clock_t endTime = clock();
	std::cout << "run time= " << (endTime - startTime) / CLOCKS_PER_SEC << 's' << endl;
}
