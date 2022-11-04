#ifndef __TFSF_H__
#define __TFSF_H__
////////////////////////////////////////////////////////////////////////////////////
//PURPOSE:
//	Total Field / Scattered Field Implementation for 3 - D FDTD code.  
//NOTES:
//	1��update1D_Einc�����ձ߽����������ַ�����������update1D_Einc()���ò�ͬ��flagѡ��ͬ�ļ��㷽����
//Reference:
//	[1]  ��±�, ����. (2011). ��Ų�ʱ�����޲�ַ����������棩.�������ӿƼ���ѧ������.
//  [2]  Taflove, A., &Hagness, S.C. (2005).Computational electrodynamics 
//       : the finite - difference time - domain method 3rd ed. Artech house.
//Record of revisions:
//   1. 06 May 2020 by Yang chao, Email: yangchaophy@foxmail.com
//   2. 02 November 2022 by Yang chao, Email: yangchaophy@foxmail.com
/////////////////////////////////////////////////////////////////////////////////////
#include<vector>
#include "matrix.h"

class TFSF
{
private:
	double dx, dy, dz;
	double alpha, thi, phi;
	double delta, dt;

	//�ܳ�����߽�
	int ItMin, ItMax, JtMin, JtMax, KtMin, KtMax;  
	int IsMin, IsMax, JsMin, JsMax, KsMin, KsMax;
	int IncidentStart, IncidentEnd, Isource;
	
	//���ӱ߽磨�ܳ��߽磩�����䳡������
	Matrix<double> Ein, Hin;
	
	std::vector<double> AttnE;
	std::vector<double> AttnH;
public:
	TFSF();
	TFSF(double dx_, double dy_, double dz_, double alpha_, double thi_, double phi_, double dt_,
		int ItMin_, int ItMax_, int JtMin_, int JtMax_, int KtMin_, int KtMax_,
		int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_,
		int IncidentStart_, int IncidentEnd_);
	
	void initializeTSFS();
	void initCoefficientsTSFS();

	void addSource(int nt);

	void update1D_Einc();
	void update1D_Hinc();
	

	void add_TFSF_Box_E(Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez);
	void add_TFSF_Box_H(Matrix<double> &Hx, Matrix<double> &Hy, Matrix<double> &Hz);
	
	void add_TFSF_X1_E(Matrix<double> &Ey, Matrix<double> &Ez);
	void add_TFSF_X2_E(Matrix<double> &Ey, Matrix<double> &Ez);
	void add_TFSF_Y1_E(Matrix<double> &Ex, Matrix<double> &Ez);
	void add_TFSF_Y2_E(Matrix<double> &Ex, Matrix<double> &Ez);
	void add_TFSF_Z1_E(Matrix<double> &Ex, Matrix<double> &Ey);
	void add_TFSF_Z2_E(Matrix<double> &Ex, Matrix<double> &Ey);


	void add_TFSF_X1_H(Matrix<double> &Hy, Matrix<double> &Hz);
	void add_TFSF_X2_H(Matrix<double> &Hy, Matrix<double> &Hz);
	void add_TFSF_Y1_H(Matrix<double> &Hx, Matrix<double> &Hz);
	void add_TFSF_Y2_H(Matrix<double> &Hx, Matrix<double> &Hz);
	void add_TFSF_Z1_H(Matrix<double> &Hx, Matrix<double> &Hy);
	void add_TFSF_Z2_H(Matrix<double> &Hx, Matrix<double> &Hy);

/////////////////////////////////////////////////////////////////
//////////////      CPML�������䲨����
////////////////////////////////////////////////////////////////
	void add_TFSF_Box_E_CPML(const Matrix<double> &ce_x_1, const Matrix<double> &ce_x_2, const Matrix<double> &ce_y_1, 
		const Matrix<double> &ce_y_2, const Matrix<double> &ce_z_1, const Matrix<double> &ce_z_2,
		Matrix<double> &psi_Exy_1, Matrix<double> &psi_Exy_2, Matrix<double> &psi_Exz_1, Matrix<double> &psi_Exz_2, Matrix<double> &psi_Eyx_1, Matrix<double> &psi_Eyx_2,
		Matrix<double> &psi_Eyz_1, Matrix<double> &psi_Eyz_2, Matrix<double> &psi_Ezx_1, Matrix<double> &psi_Ezx_2, Matrix<double> &psi_Ezy_1, Matrix<double> &psi_Ezy_2);

	void add_TFSF_Box_H_CPML(const Matrix<double> &ch_x_1, const Matrix<double> &ch_x_2, const Matrix<double> &ch_y_1, 
		const Matrix<double> &ch_y_2, const Matrix<double> &ch_z_1, const Matrix<double> &ch_z_2,
		Matrix<double> &psi_Hxy_1, Matrix<double> &psi_Hxy_2, Matrix<double> &psi_Hxz_1, Matrix<double> &psi_Hxz_2, Matrix<double> &psi_Hyx_1, Matrix<double> &psi_Hyx_2,
		Matrix<double> &psi_Hyz_1, Matrix<double> &psi_Hyz_2, Matrix<double> &psi_Hzx_1, Matrix<double> &psi_Hzx_2, Matrix<double> &psi_Hzy_1, Matrix<double> &psi_Hzy_2);

	void add_TFSF_X1_E_CPML(const Matrix<double> &ce_x_1, Matrix<double> &psi_Eyx_1, Matrix<double> &psi_Ezx_1);
	void add_TFSF_X2_E_CPML(const Matrix<double> &ce_x_2, Matrix<double> &psi_Eyx_2, Matrix<double> &psi_Ezx_2);
	void add_TFSF_Y1_E_CPML(const Matrix<double> &ce_y_1, Matrix<double> &psi_Exy_1, Matrix<double> &psi_Ezy_1);
	void add_TFSF_Y2_E_CPML(const Matrix<double> &ce_y_2, Matrix<double> &psi_Exy_2, Matrix<double> &psi_Ezy_2);
	void add_TFSF_Z1_E_CPML(const Matrix<double> &ce_z_1, Matrix<double> &psi_Exz_1, Matrix<double> &psi_Eyz_1);
	void add_TFSF_Z2_E_CPML(const Matrix<double> &ce_z_2, Matrix<double> &psi_Exz_2, Matrix<double> &psi_Eyz_2);

	void add_TFSF_X1_H_CPML(const Matrix<double> &ch_x_1, Matrix<double> &psi_Hyx_1, Matrix<double> &psi_Hzx_1);
	void add_TFSF_X2_H_CPML(const Matrix<double> &ch_x_2, Matrix<double> &psi_Hyx_2, Matrix<double> &psi_Hzx_2);
	void add_TFSF_Y1_H_CPML(const Matrix<double> &ch_y_1, Matrix<double> &psi_Hxy_1, Matrix<double> &psi_Hzy_1);
	void add_TFSF_Y2_H_CPML(const Matrix<double> &ch_y_2, Matrix<double> &psi_Hxy_2, Matrix<double> &psi_Hzy_2);
	void add_TFSF_Z1_H_CPML(const Matrix<double> &ch_z_1, Matrix<double> &psi_Hxz_1, Matrix<double> &psi_Hyz_1);
	void add_TFSF_Z2_H_CPML(const Matrix<double> &ch_z_2, Matrix<double> &psi_Hxz_2, Matrix<double> &psi_Hyz_2);

private:
	double ComputeComponentFrac();
	int source_position();
	double source(double time);
	double attenuationFactor(double x, double y, double z);
};

#endif 