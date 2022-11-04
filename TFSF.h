#pragma once
//PURPOSE:
//	Total Field / Scattered Field Implementation for 3 - D FDTD code.  
//NOTES:
//
//Reference:
//	[1]  葛德彪, 闫玉波. (2011). 电磁波时域有限差分方法（第三版）.西安电子科技大学出版社.
//  [2]  Taflove, A., &Hagness, S.C. (2005).Computational electrodynamics 
//       : the finite - difference time - domain method 3rd ed. Artech house.
//Record of revisions:
//   1. 06 May 2020 by Yang chao, Email: yangchaophy@foxmail.com
//      Northwest Institute of Nuclear Technology, Xi'an; 710024, China

#include <fstream>
#include<iostream>
#include"physical_constants.h"
#include "matrix.h"
using namespace std;
class TFSF
{
private:
	double dx, dy, dz;
	double alpha, thi, phi;
	double delta, dt;

	//总场区域边界
	int ItMin, ItMax, JtMin, JtMax, KtMin, KtMax;
	int IsMin, IsMax, JsMin, JsMax, KsMin, KsMax;
	
	int IncidentStart, IncidentEnd, Isource;

	//连接边界（总场边界）处入射场各分量
	Matrix<double> Ein, Hin;

	double AttnE[NApml];
	double AttnH[NApml];


public:
	TFSF();
	TFSF(double dx_, double dy_, double dz_, double alpha_, double thi_, double phi_, double dt_, int IncidentStart_, int IncidentEnd_,int ItMin_, 
		int ItMax_, int JtMin_, int JtMax_, int KtMin_, int KtMax_, int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_);

	void update1D_Einc(double nt);
	void update1D_Hinc();

	void add_TFSF_Box_H(Matrix<double> &den_hx, Matrix<double> &den_hy, Matrix<double> &den_hz,
		Matrix<double> &Hx, Matrix<double> &Hy, Matrix<double> &Hz);
	void add_TFSF_Box_E(double CB[][MediaNo][MediaNo][MediaNo], Matrix<int> &ob,Matrix<double> &den_ex, Matrix<double> &den_ey, 
		Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez);

	void add_TFSF_X1_H(Matrix<double> &den_hx, Matrix<double> &Hy, Matrix<double> &Hz);
	void add_TFSF_X2_H(Matrix<double> &den_hx, Matrix<double> &Hy, Matrix<double> &Hz);
	void add_TFSF_Y1_H(Matrix<double> &den_hy, Matrix<double> &Hx, Matrix<double> &Hz);
	void add_TFSF_Y2_H(Matrix<double> &den_hy, Matrix<double> &Hx, Matrix<double> &Hz);
	void add_TFSF_Z1_H(Matrix<double> &den_hz, Matrix<double> &Hx, Matrix<double> &Hy);
	void add_TFSF_Z2_H(Matrix<double> &den_hz, Matrix<double> &Hx, Matrix<double> &Hy);

	void add_TFSF_X1_E(double CB[][MediaNo][MediaNo][MediaNo],Matrix<int> &ob, Matrix<double> &den_ex, Matrix<double> &Ey, Matrix<double> &Ez);
	void add_TFSF_X2_E(double CB[][MediaNo][MediaNo][MediaNo],Matrix<int> &ob, Matrix<double> &den_ex, Matrix<double> &Ey, Matrix<double> &Ez);
	void add_TFSF_Y1_E(double CB[][MediaNo][MediaNo][MediaNo],Matrix<int> &ob, Matrix<double> &den_ey, Matrix<double> &Ex, Matrix<double> &Ez);
	void add_TFSF_Y2_E(double CB[][MediaNo][MediaNo][MediaNo],Matrix<int> &ob, Matrix<double> &den_ey, Matrix<double> &Ex, Matrix<double> &Ez);
	void add_TFSF_Z1_E(double CB[][MediaNo][MediaNo][MediaNo],Matrix<int> &ob, Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey);
	void add_TFSF_Z2_E(double CB[][MediaNo][MediaNo][MediaNo],Matrix<int> &ob, Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////      CPML区域入射波加入
	void add_TFSF_Box_E_CPML(Matrix<double> &ce_x_1, Matrix<double> &ce_x_2, Matrix<double> &ce_y_1, Matrix<double> &ce_y_2, Matrix<double> &ce_z_1, Matrix<double> &ce_z_2,
		Matrix<double> &psi_Exy_1, Matrix<double> &psi_Exy_2, Matrix<double> &psi_Exz_1, Matrix<double> &psi_Exz_2, Matrix<double> &psi_Eyx_1, Matrix<double> &psi_Eyx_2,
		Matrix<double> &psi_Eyz_1, Matrix<double> &psi_Eyz_2, Matrix<double> &psi_Ezx_1, Matrix<double> &psi_Ezx_2, Matrix<double> &psi_Ezy_1, Matrix<double> &psi_Ezy_2);

	void add_TFSF_Box_H_CPML(Matrix<double> &ch_x_1, Matrix<double> &ch_x_2, Matrix<double> &ch_y_1, Matrix<double> &ch_y_2, Matrix<double> &ch_z_1, Matrix<double> &ch_z_2,
		Matrix<double> &psi_Hxy_1, Matrix<double> &psi_Hxy_2, Matrix<double> &psi_Hxz_1, Matrix<double> &psi_Hxz_2, Matrix<double> &psi_Hyx_1, Matrix<double> &psi_Hyx_2,
		Matrix<double> &psi_Hyz_1, Matrix<double> &psi_Hyz_2, Matrix<double> &psi_Hzx_1, Matrix<double> &psi_Hzx_2, Matrix<double> &psi_Hzy_1, Matrix<double> &psi_Hzy_2);

	void add_TFSF_X1_H_CPML(Matrix<double> &ch_x_1, Matrix<double> &psi_Hyx_1, Matrix<double> &psi_Hzx_1);
	void add_TFSF_X2_H_CPML(Matrix<double> &ch_x_2, Matrix<double> &psi_Hyx_2, Matrix<double> &psi_Hzx_2);
	void add_TFSF_Y1_H_CPML(Matrix<double> &ch_y_1, Matrix<double> &psi_Hxy_1, Matrix<double> &psi_Hzy_1);
	void add_TFSF_Y2_H_CPML(Matrix<double> &ch_y_2, Matrix<double> &psi_Hxy_2, Matrix<double> &psi_Hzy_2);
	void add_TFSF_Z1_H_CPML(Matrix<double> &ch_z_1, Matrix<double> &psi_Hxz_1, Matrix<double> &psi_Hyz_1);
	void add_TFSF_Z2_H_CPML(Matrix<double> &ch_z_2, Matrix<double> &psi_Hxz_2, Matrix<double> &psi_Hyz_2);

	void add_TFSF_X1_E_CPML(Matrix<double> &ce_x_1, Matrix<double> &psi_Eyx_1, Matrix<double> &psi_Ezx_1);
	void add_TFSF_X2_E_CPML(Matrix<double> &ce_x_2, Matrix<double> &psi_Eyx_2, Matrix<double> &psi_Ezx_2);
	void add_TFSF_Y1_E_CPML(Matrix<double> &ce_y_1, Matrix<double> &psi_Exy_1, Matrix<double> &psi_Ezy_1);
	void add_TFSF_Y2_E_CPML(Matrix<double> &ce_y_2, Matrix<double> &psi_Exy_2, Matrix<double> &psi_Ezy_2);
	void add_TFSF_Z1_E_CPML(Matrix<double> &ce_z_1, Matrix<double> &psi_Exz_1, Matrix<double> &psi_Eyz_1);
	void add_TFSF_Z2_E_CPML(Matrix<double> &ce_z_2, Matrix<double> &psi_Exz_2, Matrix<double> &psi_Eyz_2);

private:
	double ComputeComponentFrac();
	int source_position();
	double source(double time);

	void APML();
};