//PURPOSE:
//	 3 - D FDTD code with CPML absorbing boundary conditions
//   This code implements the finite-difference time-domain solution of Maxwell's curl equations 
//   over a three-dimensional Cartesian space lattice. The grid is terminated by CPML absorbing
//   boundary conditions. The number of grid cells as well as the thickness of
//   the PML in each Cartesian direction can be varied independently.
//NOTES:
//	1. ����ı����FDTD.cpp�еĺ���initCoeficients���á�
//	2. �ֲ���Ĳ�����FDTD.cpp�еĺ���initCoeficients���á�
//	3. ���䲨��TSFS.cpp�еĺ���source���á�
//Reference:
//	[1]  ��±�, ����. (2011). ��Ų�ʱ�����޲�ַ����������棩.�������ӿƼ���ѧ������.
//  [2]  Taflove, A., &Hagness, S.C. (2005).Computational electrodynamics 
//       : the finite - difference time - domain method 3rd ed. Artech house.
//Record of revisions:
//   1. 06 June 2020 by Yang chao, Email: yangchaophy@foxmail.com
//      Northwest Institute of Nuclear Technology, Xi'an; 710024, China
//	 2. 19 August 2020 by Yang chao
//      �������������Ų����ĵ�Чֵ��ƽ��ֵ����ref [1]section 2.5, section 8.4.1��
#include <iostream>
#include <ctime>
#include "physical_constants.h" 
#include "FDTD.h"


int main()
{
	clock_t startTime = clock();

	//Specify Grid Cell Size in Each Direction and Calculate the
	//Resulting Courant - Stable Time Step
	double dx, dy, dz, dt;
	//cell size in each direction
	dx =0.05, dy = dx, dz = dx;   
	//time step increment
	//dt = 0.99 / (c0*sqrt(1.0 / pow(dx, 2) + 1.0 / pow(dy, 2) + 1.0 / pow(dz, 2)));    //Eq(3-2-8)	
	dt = dx / 2.0 / c0;
	
	int IsMin, IsMax, JsMin, JsMax, KsMin, KsMax;  //ɢ�䳡����߽�
	IsMin = -20, IsMax = 20;
	JsMin = IsMin, JsMax = IsMax;
	KsMin =-20, KsMax = 20;
	int tmax; //total number of time steps
	tmax =300;   //<500

	//Specify the CPML Thickness in Each Direction(Value of Zero
	//Corresponds to No PML, and the Grid is Terminated with a PEC)
	//PML thickness in each direction
	int nxPML_1, nxPML_2, nyPML_1, nyPML_2, nzPML_1, nzPML_2;
	nxPML_1 = 16, nxPML_2 = nxPML_1, nyPML_1 = nxPML_1, nyPML_2 = nxPML_1, nzPML_1 = nxPML_1, nzPML_2 = nxPML_1;
	
	//Specify Material Relative Permittivity and Conductivity
	double epsr[MediaNo] = { 1.0, 1.0 };
	double sig[MediaNo] = { 0.0, 0.0 };

	//************************************************************************************
	//ƽ�沨����
	//************************************************************************************
	int TFSFGrid = -5;
	int ItMin = IsMin - TFSFGrid;   //connective boundary
	int ItMax = IsMax + TFSFGrid;
	int JtMin = JsMin - TFSFGrid;
	int JtMax = JsMax + TFSFGrid;
	int KtMin = KsMin + 5;
	int KtMax = KsMax + 10;

	double alpha = 0.0 / 180.0*pi;
	double thi = 180.0 / 180.0*pi;
	double phi = 180.0 / 180.0*pi;
	int IncStart = -500;
	int IncEnd = 500;
	//************************************************************************************

	FDTD fdtd(tmax, IsMin, IsMax, JsMin, JsMax, KsMin, KsMax, nxPML_1, nxPML_2, nyPML_1, nyPML_2, nzPML_1, nzPML_2,
		dt, dx, dy, dz, epsr, sig, alpha, thi, phi, IncStart, IncEnd, ItMin, ItMax, JtMin, JtMax, KtMin, KtMax);
	
	fdtd.initCoeficients();
	fdtd.compute();

	clock_t endTime = clock();
	cout << "run time= " << (endTime - startTime) / CLOCKS_PER_SEC << 's' << endl;
	return 0;
}