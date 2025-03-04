#ifndef __TFSF_H__
#define __TFSF_H__
////////////////////////////////////////////////////////////////////////////////////
//PURPOSE:
//	Total Field / Scattered Field Implementation for 3 - D FDTD code.  
//NOTES:
//	1、update1D_Einc中吸收边界条件有两种方法，可以在update1D_Einc()设置不同的flag选择不同的计算方法。
//Reference:
//	[1]  葛德彪, 闫玉波. (2011). 电磁波时域有限差分方法（第三版）.西安电子科技大学出版社.
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

	//总场区域边界
	int ItMin, ItMax, JtMin, JtMax, KtMin, KtMax;  
	int IsMin, IsMax, JsMin, JsMax, KsMin, KsMax;
	int IncidentStart, IncidentEnd, Isource;
	
	//连接边界（总场边界）处入射场各分量
	Matrix<double> Ein, Hin;
	
	std::vector<double> AttnE00;
	std::vector<double> AttnH00;
	std::vector<double> AttnE30;
	std::vector<double> AttnH30;
	std::vector<double> AttnE45;
	std::vector<double> AttnH45;
	std::vector<double> AttnE60;
	std::vector<double> AttnH60;
	std::vector<double> AttnE90;
	std::vector<double> AttnH90;

	std::vector<double> AttnEx;
	std::vector<double> AttnHx;
	std::vector<double> AttnEy;
	std::vector<double> AttnHy;
	std::vector<double> AttnEz;
	std::vector<double> AttnHz;

	std::vector<double> phaseE00;
	std::vector<double> phaseH00;
	std::vector<double> phaseE30;
	std::vector<double> phaseH30;
	std::vector<double> phaseE45;
	std::vector<double> phaseH45;
	std::vector<double> phaseE60;
	std::vector<double> phaseH60;
	std::vector<double> phaseE90;
	std::vector<double> phaseH90;

	std::vector<double> phaseEx;
	std::vector<double> phaseHx;
	std::vector<double> phaseEy;
	std::vector<double> phaseHy;
	std::vector<double> phaseEz;
	std::vector<double> phaseHz;
	
	////bAnalysis=true解析计算，bAnalysis=false数值计算。
	////应该选择bAnalysis=false数值计算。
	bool bAnalysis;
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
	

	void add_TFSF_Box_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ex, const Matrix<double> &den_ey,
		const Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez);
	void add_TFSF_Box_H(const Matrix<double> &den_hx, const Matrix<double> &den_hy, const Matrix<double> &den_hz,
		Matrix<double> &Hx, Matrix<double> &Hy, Matrix<double> &Hz);
	
	void add_TFSF_X1_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ex, Matrix<double> &Ey, Matrix<double> &Ez);
	void add_TFSF_X2_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ex, Matrix<double> &Ey, Matrix<double> &Ez);
	void add_TFSF_Y1_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ey, Matrix<double> &Ex, Matrix<double> &Ez);
	void add_TFSF_Y2_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ey, Matrix<double> &Ex, Matrix<double> &Ez);
	void add_TFSF_Z1_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey);
	void add_TFSF_Z2_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey);


	void add_TFSF_X1_H(const Matrix<double> &den_hx, Matrix<double> &Hy, Matrix<double> &Hz);
	void add_TFSF_X2_H(const Matrix<double> &den_hx, Matrix<double> &Hy, Matrix<double> &Hz);
	void add_TFSF_Y1_H(const Matrix<double> &den_hy, Matrix<double> &Hx, Matrix<double> &Hz);
	void add_TFSF_Y2_H(const Matrix<double> &den_hy, Matrix<double> &Hx, Matrix<double> &Hz);
	void add_TFSF_Z1_H(const Matrix<double> &den_hz, Matrix<double> &Hx, Matrix<double> &Hy);
	void add_TFSF_Z2_H(const Matrix<double> &den_hz, Matrix<double> &Hx, Matrix<double> &Hy);

/////////////////////////////////////////////////////////////////
//////////////      CPML区域入射波加入
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

	/////////////////////////////////////////////////////////////////
	//////////////    解析法
	////////////////////////////////////////////////////////////////
	void add_TSFS_Box_E_analysis(double nt, const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ex,
		const Matrix<double> &den_ey, const Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez);
	void add_TSFS_Box_H_analysis(double nt, const Matrix<double> &den_hx, const Matrix<double> &den_hy, const Matrix<double> &den_hz, 
		Matrix<double> &Hx, Matrix<double> &Hy, Matrix<double> &Hz);


	void add_TFSF_Box_E_CPML_analysis(double nt, const Matrix<double> &ce_x_1, const Matrix<double> &ce_x_2, const Matrix<double> &ce_y_1,
		const Matrix<double> &ce_y_2, const Matrix<double> &ce_z_1, const Matrix<double> &ce_z_2,
		Matrix<double> &psi_Exy_1, Matrix<double> &psi_Exy_2, Matrix<double> &psi_Exz_1, Matrix<double> &psi_Exz_2, Matrix<double> &psi_Eyx_1, Matrix<double> &psi_Eyx_2,
		Matrix<double> &psi_Eyz_1, Matrix<double> &psi_Eyz_2, Matrix<double> &psi_Ezx_1, Matrix<double> &psi_Ezx_2, Matrix<double> &psi_Ezy_1, Matrix<double> &psi_Ezy_2);
	
	void add_TFSF_Box_H_CPML_analysis(double nt, const Matrix<double> &ch_x_1, const Matrix<double> &ch_x_2, const Matrix<double> &ch_y_1,
		const Matrix<double> &ch_y_2, const Matrix<double> &ch_z_1, const Matrix<double> &ch_z_2,
		Matrix<double> &psi_Hxy_1, Matrix<double> &psi_Hxy_2, Matrix<double> &psi_Hxz_1, Matrix<double> &psi_Hxz_2, Matrix<double> &psi_Hyx_1, Matrix<double> &psi_Hyx_2,
		Matrix<double> &psi_Hyz_1, Matrix<double> &psi_Hyz_2, Matrix<double> &psi_Hzx_1, Matrix<double> &psi_Hzx_2, Matrix<double> &psi_Hzy_1, Matrix<double> &psi_Hzy_2);

	void add_TFSF_X1_E_CPML_analysis(double nt, const Matrix<double> &ce_x_1, Matrix<double> &psi_Eyx_1, Matrix<double> &psi_Ezx_1);
	void add_TFSF_X2_E_CPML_analysis(double nt, const Matrix<double> &ce_x_2, Matrix<double> &psi_Eyx_2, Matrix<double> &psi_Ezx_2);
	void add_TFSF_Y1_E_CPML_analysis(double nt, const Matrix<double> &ce_y_1, Matrix<double> &psi_Exy_1, Matrix<double> &psi_Ezy_1);
	void add_TFSF_Y2_E_CPML_analysis(double nt, const Matrix<double> &ce_y_2, Matrix<double> &psi_Exy_2, Matrix<double> &psi_Ezy_2);
	void add_TFSF_Z1_E_CPML_analysis(double nt, const Matrix<double> &ce_z_1, Matrix<double> &psi_Exz_1, Matrix<double> &psi_Eyz_1);
	void add_TFSF_Z2_E_CPML_analysis(double nt, const Matrix<double> &ce_z_2, Matrix<double> &psi_Exz_2, Matrix<double> &psi_Eyz_2);

	void add_TFSF_X1_H_CPML_analysis(double nt, const Matrix<double> &ch_x_1, Matrix<double> &psi_Hyx_1, Matrix<double> &psi_Hzx_1);
	void add_TFSF_X2_H_CPML_analysis(double nt, const Matrix<double> &ch_x_2, Matrix<double> &psi_Hyx_2, Matrix<double> &psi_Hzx_2);
	void add_TFSF_Y1_H_CPML_analysis(double nt, const Matrix<double> &ch_y_1, Matrix<double> &psi_Hxy_1, Matrix<double> &psi_Hzy_1);
	void add_TFSF_Y2_H_CPML_analysis(double nt, const Matrix<double> &ch_y_2, Matrix<double> &psi_Hxy_2, Matrix<double> &psi_Hzy_2);
	void add_TFSF_Z1_H_CPML_analysis(double nt, const Matrix<double> &ch_z_1, Matrix<double> &psi_Hxz_1, Matrix<double> &psi_Hyz_1);
	void add_TFSF_Z2_H_CPML_analysis(double nt, const Matrix<double> &ch_z_2, Matrix<double> &psi_Hxz_2, Matrix<double> &psi_Hyz_2);



private:
	double ComputeComponentFrac();
	int source_position();
	double source(double time, double phase=0.0);
	double attenuationFactor(double x, double y, double z);
	double phaseCorrection(double x, double y, double z);
	void CalAttnEH();
	void CalPhaseEH();
};

#endif 