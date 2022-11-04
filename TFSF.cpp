#include<iostream>
#include <fstream>
#include<string>
#include<cmath>
#include"physical_constants.h"
#include"TFSF.h"


TFSF::TFSF()
{
	std::cout << "TSFS: void constructor\n";
}


TFSF::TFSF(double dx_, double dy_, double dz_, double alpha_, double thi_, double phi_, double dt_, 
	int ItMin_, int ItMax_, int JtMin_, int JtMax_, int KtMin_, int KtMax_,
	int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_,
	int IncidentStart_,int IncidentEnd_)
{
	dx = dx_;
	dy = dy_;
	dz = dz_;

	alpha = alpha_;  //0=垂直极化，90=水平极化
	thi = thi_;
	phi = phi_;

	delta = ComputeComponentFrac()*dx;
	dt = dt_;

	ItMin = ItMin_;
	ItMax = ItMax_;
	JtMin = JtMin_;
	JtMax = JtMax_;
	KtMin = KtMin_;
	KtMax = KtMax_;

	IsMin = IsMin_;
	IsMax = IsMax_;
	JsMin = JsMin_;
	JsMax = JsMax_;
	KsMin = KsMin_;
	KsMax = KsMax_;

	IncidentStart = IncidentStart_;  // 平行x轴， j=Jtmin
	IncidentEnd = IncidentEnd_;
	Isource = source_position();
}

/////////////////////////////////////////////////////////////////////////////////////
//PURPOSE:
//	计算一维辅助入射波空间离散间隔与三维离散间隔的比值：Delta/dx
//Reference:
//	[1]  葛德彪, 闫玉波. (2011). 电磁波时域有限差分方法（第三版）.西安电子科技大学出版社.
//  [2]  Taflove, A., &Hagness, S.C. (2005).Computational electrodynamics 
//       : the finite - difference time - domain method 3rd ed. Artech house.
//////////////////////////////////////////////////////////////////////////////////////
double TFSF::ComputeComponentFrac()
{
	double st = sin(thi);
	double ct = cos(thi);
	double sp = sin(phi);
	double cp = cos(phi);
	
	//ref[2], P217, Eq(5.69)
	double m_VFrac = sqrt(pow(st, 4)*(pow(cp, 4) + pow(sp, 4)) + pow(ct, 4));
	
	//std::cout << "m_VFrac= " << m_VFrac << std::endl;
	return m_VFrac;
}

int TFSF::source_position()
{
	int i0, j0, k0;
	if (thi <= pi / 2.0) {
		if (phi <= pi / 2.0) {
			i0 = ItMin;
			j0 = JtMin;
			k0 = KtMin;
		}
		else if (phi <= pi) {
			i0 = ItMax;
			j0 = JtMin;
			k0 = KtMin;
		}
		else if (phi <= 1.5*pi) {
			i0 = ItMax;
			j0 = JtMax;
			k0 = KtMin;
		}
		else {
			i0 = ItMin;
			j0 = JtMax;
			k0 = KtMin;
		}
	}
	else {
		if (phi <= pi / 2.0) {
			i0 = ItMin;
			j0 = JtMin;
			k0 = KtMax;
		}
		else if (phi <= pi) {
			i0 = ItMax;
			j0 = JtMin;
			k0 = KtMax;
		}
		else if (phi <= 1.5*pi) {
			i0 = ItMax;
			j0 = JtMax;
			k0 = KtMax;
		}
		else {
			i0 = ItMin;
			j0 = JtMax;
			k0 = KtMax;
		}
	}

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);

	int Is = int(i0*kx + j0*ky + k0*kz) - 2;
	return Is;
}


void TFSF::initializeTSFS()
{
	Ein.Allocate(IncidentStart, IncidentEnd);
	Hin.Allocate(IncidentStart, IncidentEnd - 1);
}


void TFSF::initCoefficientsTSFS()
{
	std::string fname = "../read_data/Z_Gauss300MHz_dx=2.5cm_m=4_a=0_k=1_N=16.dat";
	std::ifstream fin(fname);

	if (!fin.is_open()) {
		std::cout << fname << ":  open fail" << std::endl;
	}

	int temp;
	double attnEtemp, attnHtemp;
	while (fin >> temp && fin >> attnEtemp&& fin >> attnHtemp)
	{
		AttnE.push_back(attnEtemp);
		AttnH.push_back(attnHtemp);
	}

	fin.close();
}

void TFSF::update1D_Einc()
{
	double FE = dt / eps0 / delta;  //一维无耗介质  sig=0，sig_m=0
	double EBin[4];
	//上一时刻的边界点。
	EBin[0] = Ein(IncidentStart);  //上一时刻
	EBin[1] = Ein(IncidentStart + 1); //
	EBin[2] = Ein(IncidentEnd - 1);
	EBin[3] = Ein(IncidentEnd);

	for (int i = IncidentStart + 1; i <= IncidentEnd - 1; i++) {
		Ein(i) = Ein(i) - FE*(Hin(i) - Hin(i - 1));
	}
	int flag = 0;
	if (flag == 0) {//一维入射波吸收边界 Eq(4-3-17)
		Ein(IncidentStart) = EBin[1] + (c0*dt - delta) / (c0*dt + delta)*(Ein(IncidentStart + 1) - EBin[0]);    
		Ein(IncidentEnd) = EBin[2] + (c0*dt - delta) / (c0*dt + delta)*(Ein(IncidentEnd - 1) - EBin[3]);
	}
	else {//一维入射波吸收边界 Eq(4-3-20)
		Ein(IncidentStart) = EBin[0] - c0*dt / delta*(EBin[0] - EBin[1]);    
		Ein(IncidentEnd) = EBin[3] - c0*dt / delta*(EBin[3] - EBin[2]);
	}
}

void TFSF::update1D_Hinc()
{
	double FH = dt / mu0 / delta;
	// 产生入射波  Ein  Hin
	for (int i = IncidentStart; i <= IncidentEnd - 1; i++) {
		Hin(i) = Hin(i) - FH*(Ein(i + 1) - Ein(i));
	}
}


void TFSF::add_TFSF_Box_E(Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez)
{
	add_TFSF_X1_E(Ey, Ez);
	add_TFSF_X2_E(Ey, Ez);
	add_TFSF_Y1_E(Ex, Ez);
	add_TFSF_Y2_E(Ex, Ez);
	add_TFSF_Z1_E(Ex, Ey);
	add_TFSF_Z2_E(Ex, Ey);
}

void TFSF::add_TFSF_X1_E(Matrix<double> &Ey, Matrix<double> &Ez)
{
	double APML = 1.0;

	int j, k;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-18)
	double Hytemp = cos(phi)*cos(alpha) - cos(thi)*sin(phi)*sin(alpha);
	double Hztemp = sin(thi)*sin(alpha);

	//************************************************************************
	//***************      x1        ***************
	//    x1  Ey
	double HzCB_x1;
	for (j = JtMin; j <= JtMax - 1; j++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = (ItMin - 0.5)*kx + (j + 0.5)*ky + k*kz;  //Hy
			II = int(T1 - 0.5);
			if (T1 >= 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HzCB_x1 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HzCB_x1 *= Hztemp;

			APML= attenuationFactor(ItMin - 0.5, j + 0.5, k);
			Ey(ItMin, j, k) += dt / eps0 / dx*HzCB_x1*APML;   //表6-3 修改mu0
		}
	}
	//   x1  Ez
	double HyCB_x1;
	for (j = JtMin; j <= JtMax; j++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = (ItMin - 0.5)*kx + j *ky + (k + 0.5)*kz;   //Hx
			II = int(T1 - 0.5);
			if (T1 > 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HyCB_x1 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HyCB_x1 *= Hytemp;

			APML = attenuationFactor(ItMin - 0.5, j, k + 0.5);
			Ez(ItMin, j, k) -= dt / eps0 / dx*HyCB_x1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_X2_E(Matrix<double> &Ey, Matrix<double> &Ez)
{
	double APML = 1.0;

	int j, k;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-18)
	double Hytemp = cos(phi)*cos(alpha) - cos(thi)*sin(phi)*sin(alpha);
	double Hztemp = sin(thi)*sin(alpha);

	//************************************************************************
	//***************      x2        ***************
	//    x2  Ey
	double HzCB_x2;
	for (j = JtMin; j <= JtMax - 1; j++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = (ItMax + 0.5)*kx + (j + 0.5)*ky + k*kz;  //Hy
			II = int(T1 - 0.5);
			if (T1 >= 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HzCB_x2 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HzCB_x2 *= Hztemp;

			APML = attenuationFactor(ItMax + 0.5, j + 0.5, k);
			Ey(ItMax, j, k) -= dt / eps0 / dx*HzCB_x2*APML;   //表6-3 修改mu0
		}
	}
	//   x2  Ez
	double HyCB_x2;
	for (j = JtMin; j <= JtMax; j++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = (ItMax + 0.5)*kx + j *ky + (k + 0.5)*kz;   //Hx
			II = int(T1 - 0.5);
			if (T1 > 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HyCB_x2 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HyCB_x2 *= Hytemp;

			APML = attenuationFactor(ItMax + 0.5, j, k + 0.5);
			Ez(ItMax, j, k) += dt / eps0 / dx*HyCB_x2*APML;   //表6-3 修改mu0
		}
	}

}

void TFSF::add_TFSF_Y1_E(Matrix<double> &Ex, Matrix<double> &Ez)
{
	double APML = 1.0;

	int i, k;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-18)
	double Hxtemp = -sin(phi)*cos(alpha) - cos(thi)*cos(phi)*sin(alpha);
	double Hztemp = sin(thi)*sin(alpha);

	//************************************************************************
	//***************      y1        ***************
	//    y1  Ez
	double HxCB_y1;
	for (i = ItMin; i <= ItMax; i++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = i*kx + (JtMin - 0.5)*ky + (k + 0.5)*kz;  //Hx
			II = int(T1 - 0.5);
			if (T1 >= 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HxCB_y1 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HxCB_y1 *= Hxtemp;

			APML = attenuationFactor(i, JtMin - 0.5, k + 0.5);
			Ez(i, JtMin, k) += dt / eps0 / dy*HxCB_y1*APML;   //表6-3 修改mu0
		}
	}
	//    y1  Ex
	double HzCB_y1;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = (i + 0.5)*kx + (JtMin - 0.5)*ky + k *kz;  //Hz
			II = int(T1 - 0.5);
			if (T1 >= 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HzCB_y1 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HzCB_y1 *= Hztemp;

			APML = attenuationFactor(i + 0.5, JtMin - 0.5, k);
			Ex(i, JtMin, k) -= dt / eps0 / dy*HzCB_y1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y2_E(Matrix<double> &Ex, Matrix<double> &Ez)
{
	double APML = 1.0;

	int i, k;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-18)
	double Hxtemp = -sin(phi)*cos(alpha) - cos(thi)*cos(phi)*sin(alpha);
	double Hztemp = sin(thi)*sin(alpha);
	//************************************************************************
	//***************      y2        ***************
	//    y2  Ez
	double HxCB_y2;
	for (i = ItMin; i <= ItMax; i++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = i*kx + (JtMax + 0.5)*ky + (k + 0.5)*kz;  //Hx
			II = int(T1 - 0.5);
			if (T1 >= 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HxCB_y2 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HxCB_y2 *= Hxtemp;

			APML = attenuationFactor(i, JtMax + 0.5, k + 0.5);
			Ez(i, JtMax, k) -= dt / eps0 / dy*HxCB_y2*APML;   //表6-3 修改mu0
		}
	}
	//    y2  Ex
	double HzCB_y2;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = (i + 0.5)*kx + (JtMax + 0.5)*ky + k *kz;  //Hz
			II = int(T1 - 0.5);
			if (T1 >= 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HzCB_y2 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HzCB_y2 *= Hztemp;

			APML = attenuationFactor(i + 0.5, JtMax + 0.5, k);
			Ex(i, JtMax, k) += dt / eps0 / dy*HzCB_y2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z1_E(Matrix<double> &Ex, Matrix<double> &Ey)
{
	double APML = 1.0;

	int i, j;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-18)
	double Hxtemp = -sin(phi)*cos(alpha) - cos(thi)*cos(phi)*sin(alpha);
	double Hytemp = cos(phi)*cos(alpha) - cos(thi)*sin(phi)*sin(alpha);
	//************************************************************************
	//***************      z1        ***************
	//    z1  Ex
	double HyCB_z1;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (j = JtMin; j <= JtMax; j++) {
			T1 = (i + 0.5)*kx + j*ky + (KtMin - 0.5)*kz;  //Hy
			II = int(T1 - 0.5);
			if (T1 >= 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HyCB_z1 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HyCB_z1 *= Hytemp;

			APML = attenuationFactor(i + 0.5, j, KtMin - 0.5);
			Ex(i, j, KtMin) += dt / eps0 / dz*HyCB_z1*APML;   //表6-3 修改mu0
		}
	}
	//   z2  Ey
	double HxCB_z1;
	for (i = ItMin; i <= ItMax; i++) {
		for (j = JtMin; j <= JtMax - 1; j++) {
			T1 = i*kx + (j + 0.5)*ky + (KtMin - 0.5)*kz;   //Hx
			II = int(T1 - 0.5);
			if (T1 > 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HxCB_z1 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HxCB_z1 *= Hxtemp;

			APML = attenuationFactor(i, j + 0.5, KtMin - 0.5);
			Ey(i, j, KtMin) -= dt / eps0 / dz*HxCB_z1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z2_E(Matrix<double> &Ex, Matrix<double> &Ey)
{
	double APML = 1.0;

	int i, j;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-18)
	double Hxtemp = -sin(phi)*cos(alpha) - cos(thi)*cos(phi)*sin(alpha);
	double Hytemp = cos(phi)*cos(alpha) - cos(thi)*sin(phi)*sin(alpha);

	//************************************************************************
	//***************      z2        ***************
	//    z2  Ex
	double HyCB_z2;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (j = JtMin; j <= JtMax; j++) {
			T1 = (i + 0.5)*kx + j*ky + (KtMax + 0.5)*kz;  //Hy
			II = int(T1 - 0.5);
			if (T1 >= 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HyCB_z2 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HyCB_z2 *= Hytemp;

			APML = attenuationFactor(i + 0.5, j, KtMax + 0.5);
			Ex(i, j, KtMax) -= dt / eps0 / dz*HyCB_z2*APML;   //表6-3 修改mu0
		}
	}
	//   z2  Ey
	double HxCB_z2;
	for (i = ItMin; i <= ItMax; i++) {
		for (j = JtMin; j <= JtMax - 1; j++) {
			T1 = i*kx + (j + 0.5)*ky + (KtMax + 0.5)*kz;   //Hx
			II = int(T1 - 0.5);
			if (T1 > 0.5) {
				III = II + 1;
				T1 = double(III) + 0.5 - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - 0.5 - double(III);
			}
			HxCB_z2 = T1*(Hin(II) - Hin(III)) + Hin(III);
			HxCB_z2 *= Hxtemp;

			APML = attenuationFactor(i, j + 0.5, KtMax + 0.5);
			Ey(i, j, KtMax) += dt / eps0 / dz*HxCB_z2*APML;   //表6-3 修改mu0
		}
	}
}


void TFSF::add_TFSF_Box_H(Matrix<double> &Hx, Matrix<double> &Hy, Matrix<double> &Hz)
{
	add_TFSF_X1_H(Hy, Hz);
	add_TFSF_X2_H(Hy, Hz);
	add_TFSF_Y1_H(Hx, Hz);
	add_TFSF_Y2_H(Hx, Hz);
	add_TFSF_Z1_H(Hx, Hy);
	add_TFSF_Z2_H(Hx, Hy);
}

void TFSF::add_TFSF_X1_H(Matrix<double> &Hy, Matrix<double> &Hz)
{
	double APML = 1.0;

	int j, k;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-17)
	double Eytemp = cos(phi)*sin(alpha) + cos(thi)*sin(phi)*cos(alpha);
	double Eztemp = -sin(thi)*cos(alpha);

	//****************************************************************************************
	//***************       x1      ***************
	//x1 Hy
	double EzCB_x1;
	for (j = JtMin; j <= JtMax; j++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = ItMin*kx + j*ky + (k + 0.5)*kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			EzCB_x1 = T1*(Ein(II) - Ein(III)) + Ein(III);
			EzCB_x1 *= Eztemp;

			APML = attenuationFactor(ItMin, j, k + 0.5);
			Hy(ItMin - 1, j, k) -=  dt / mu0 / dx*EzCB_x1*APML;   //表6-3 修改mu0
		}
	}
	//x1 Hz
	double EyCB_x1;
	for (j = JtMin; j <= JtMax - 1; j++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = ItMin*kx + (j + 0.5)*ky + k *kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			EyCB_x1 = T1*(Ein(II) - Ein(III)) + Ein(III);
			EyCB_x1 *= Eytemp;

			APML = attenuationFactor(ItMin, j + 0.5 , k);
			Hz(ItMin - 1, j, k) += dt / mu0 / dx*EyCB_x1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_X2_H(Matrix<double> &Hy, Matrix<double> &Hz)
{
	double APML = 1.0;

	int j, k;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-17)
	double Eytemp = cos(phi)*sin(alpha) + cos(thi)*sin(phi)*cos(alpha);
	double Eztemp = -sin(thi)*cos(alpha);

	//****************************************************************************************
	//***************       x2      ***************
	//x2 Hy
	double EzCB_x2;
	for (j = JtMin; j <= JtMax; j++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = ItMax*kx + j*ky + (k + 0.5)*kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			EzCB_x2 = T1*(Ein(II) - Ein(III)) + Ein(III);
			EzCB_x2 *= Eztemp;

			APML = attenuationFactor(ItMax, j, k + 0.5);
			Hy(ItMax, j, k) +=  dt / mu0 / dx*EzCB_x2*APML;   //表6-3 修改mu0
		}
	}

	//x2 Hz
	double EyCB_x2;
	for (j = JtMin; j <= JtMax - 1; j++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = ItMax*kx + (j + 0.5)*ky + k *kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			EyCB_x2 = T1*(Ein(II) - Ein(III)) + Ein(III);
			EyCB_x2 *= Eytemp;

			APML = attenuationFactor(ItMax, j + 0.5, k);
			Hz(ItMax, j, k) -=  dt / mu0 / dx*EyCB_x2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y1_H(Matrix<double> &Hx, Matrix<double> &Hz)
{
	double APML = 1.0;

	int i, k;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-17)
	double Extemp = -sin(phi)*sin(alpha) + cos(thi)*cos(phi)*cos(alpha);
	double Eztemp = -sin(thi)*cos(alpha);
	//****************************************************************************************
	//***************       y1      ***************
	// y1 Hz
	double ExCB_y1;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = (i + 0.5)*kx + JtMin*ky + k*kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			ExCB_y1 = T1*(Ein(II) - Ein(III)) + Ein(III);
			ExCB_y1 *= Extemp;

			APML = attenuationFactor(i + 0.5, JtMin, k);
			Hz(i, JtMin - 1, k) -=  dt / mu0 / dy*ExCB_y1*APML;   //表6-3 修改mu0
		}
	}
	// y1  Hx
	double EzCB_y1;
	for (i = ItMin; i <= ItMax; i++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = i*kx + JtMin*ky + (k + 0.5)*kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			EzCB_y1 = T1*(Ein(II) - Ein(III)) + Ein(III);
			EzCB_y1 *= Eztemp;

			APML = attenuationFactor(i, JtMin, k + 0.5);
			Hx(i, JtMin - 1, k) +=  dt / mu0 / dy*EzCB_y1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y2_H(Matrix<double> &Hx, Matrix<double> &Hz)
{
	double APML = 1.0;

	int i, k;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-17)
	double Extemp = -sin(phi)*sin(alpha) + cos(thi)*cos(phi)*cos(alpha);
	double Eztemp = -sin(thi)*cos(alpha);
	//****************************************************************************************
	//***************       y2      ***************
	// y2  Hz
	double ExCB_y2;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = (i + 0.5)*kx + JtMax*ky + k*kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			ExCB_y2 = T1*(Ein(II) - Ein(III)) + Ein(III);
			ExCB_y2 *= Extemp;

			APML = attenuationFactor(i + 0.5, JtMax, k);
			Hz(i, JtMax, k) += dt / mu0 / dy*ExCB_y2*APML;   //表6-3 修改mu0
		}
	}
	// y2  Hx
	double EzCB_y2;
	for (i = ItMin; i <= ItMax; i++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = i*kx + JtMax*ky + (k + 0.5)*kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			EzCB_y2 = T1*(Ein(II) - Ein(III)) + Ein(III);
			EzCB_y2 *= Eztemp;

			APML = attenuationFactor(i, JtMax, k + 0.5);
			Hx(i, JtMax, k) -=  dt / mu0 / dy*EzCB_y2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z1_H(Matrix<double> &Hx, Matrix<double> &Hy)
{
	double APML = 1.0;

	int i, j;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-17)
	double Extemp = -sin(phi)*sin(alpha) + cos(thi)*cos(phi)*cos(alpha);
	double Eytemp = cos(phi)*sin(alpha) + cos(thi)*sin(phi)*cos(alpha);
	//****************************************************************************************
	//***************       z1      ***************
	//z1  Hx
	double EyCB_z1;
	for (i = ItMin; i <= ItMax; i++) {
		for (j = JtMin; j <= JtMax - 1; j++) {
			T1 = i*kx + (j + 0.5)*ky + KtMin*kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			EyCB_z1 = T1*(Ein(II) - Ein(III)) + Ein(III);
			EyCB_z1 *= Eytemp;

			APML = attenuationFactor(i, j + 0.5, KtMin);
			Hx(i, j, KtMin - 1) -= dt / mu0 / dz*EyCB_z1*APML;   //表6-3 修改mu0
		}
	}

	//z1  Hy
	double ExCB_z1;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (j = JtMin; j <= JtMax; j++) {
			T1 = (i + 0.5)*kx + j*ky + KtMin*kz;
			II = int(T1);  //舍弃尾部
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			ExCB_z1 = T1*(Ein(II) - Ein(III)) + Ein(III);
			ExCB_z1 *= Extemp;

			APML = attenuationFactor(i + 0.5, j, KtMin);
			Hy(i, j, KtMin - 1) += dt / mu0 / dz*ExCB_z1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z2_H(Matrix<double> &Hx, Matrix<double> &Hy)
{
	double APML = 1.0;

	int i, j;
	int II, III;
	double T1;

	double kx = dx / delta*sin(thi)*cos(phi);
	double ky = dy / delta*sin(thi)*sin(phi);
	double kz = dz / delta*cos(thi);
	//[1] Eq(6-7-17)
	double Extemp = -sin(phi)*sin(alpha) + cos(thi)*cos(phi)*cos(alpha);
	double Eytemp = cos(phi)*sin(alpha) + cos(thi)*sin(phi)*cos(alpha);

	//****************************************************************************************
	//***************       z2     ***************
	//z2  Hx
	double EyCB_z2;
	for (i = ItMin; i <= ItMax; i++) {
		for (j = JtMin; j <= JtMax - 1; j++) {
			T1 = i*kx + (j + 0.5)*ky + KtMax*kz;
			II = int(T1);
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			EyCB_z2 = T1*(Ein(II) - Ein(III)) + Ein(III);
			EyCB_z2 *= Eytemp;

			APML = attenuationFactor(i, j + 0.5, KtMax);
			Hx(i, j, KtMax) += dt / mu0 / dz*EyCB_z2*APML;   //表6-3 修改mu0
		}
	}

	//z2  Hy
	double ExCB_z2;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (j = JtMin; j <= JtMax; j++) {
			T1 = (i + 0.5)*kx + j*ky + KtMax*kz;
			II = int(T1);  //舍弃尾部
			if (T1 >= 0) {
				III = II + 1;
				T1 = double(III) - T1;
			}
			else {
				III = II - 1;
				T1 = T1 - double(III);
			}
			ExCB_z2 = T1*(Ein(II) - Ein(III)) + Ein(III);
			ExCB_z2 *= Extemp;

			APML = attenuationFactor(i + 0.5, j, KtMax);
			Hy(i, j, KtMax) -= dt / mu0 / dz*ExCB_z2*APML;   //表6-3 修改mu0
		}
	}
}


double TFSF::attenuationFactor(double x, double y, double z)
{
	int i, j, k;
	bool bIsAddHalfX, bIsAddHalfY, bIsAddHalfZ;

	////////////  i位置判断
	int i0 = ((x >= 0) ? int(x) : int(x) - 1);
	if (x - i0 < 0.25) {
		i = i0;
		bIsAddHalfX = false;
	}
	else if (x - i0 > 0.75) {
		i = i0 + 1;
		bIsAddHalfX = false;
	}
	else {
		i = i0;
		bIsAddHalfX = true;
	}

	////////////  j位置判断
	int j0 = ((y >= 0) ? int(y) : int(y) - 1);
	if (y - j0 < 0.25) {
		j = j0;
		bIsAddHalfY = false;
	}
	else if (y - j0 > 0.75) {
		j = j0 + 1;
		bIsAddHalfY = false;
	}
	else {
		j = j0;
		bIsAddHalfY = true;
	}

	////////////  k位置判断
	int k0 = ((z >= 0) ? int(z) : int(z) - 1);
	if (z - k0 < 0.25) {
		k = k0;
		bIsAddHalfZ = false;
	}
	else if (z - k0 > 0.75) {
		k = k0 + 1;
		bIsAddHalfZ = false;
	}
	else {
		k = k0;
		bIsAddHalfZ = true;
	}


	double AxPML, AyPML, AzPML, APML;
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	////////////////////////////////////////////
	//******** x 方向
	if (!bIsAddHalfX)
	{
		if (i < IsMin)
			AxPML = pow(AttnE[IsMin - i], -k0x);
		else if (i > IsMax)
			AxPML = pow(AttnE[i - IsMax], k0x);
		else
			AxPML = 1.0;
	}
	else {
		if (i < IsMin)
			AxPML = pow(AttnH[IsMin - i - 1], -k0x);
		else if (i >= IsMax)
			AxPML = pow(AttnH[i - IsMax], k0x);
		else
			AxPML = 1.0;
	}


	////////////////////////////////////////////
	//******** y 方向
	if (!bIsAddHalfY)
	{
		if (j < JsMin)
			AyPML = pow(AttnE[JsMin - j], -k0y);
		else if (j > JsMax)
			AyPML = pow(AttnE[j - JsMax], k0y);
		else
			AyPML = 1.0;
	}
	else
	{
		if (j < JsMin)
			AyPML = pow(AttnH[JsMin - j - 1], -k0y);
		else if (j >= JsMax)
			AyPML = pow(AttnH[j - JsMax], k0y);
		else
			AyPML = 1.0;
	}


	////////////////////////////////////////////
	//******** z 方向
	if (!bIsAddHalfZ)
	{
		if (k < KsMin)
			AzPML = pow(AttnE[KsMin - k], -k0z);
		else if (k > KsMax)
			AzPML = pow(AttnE[k - KsMax], k0z);
		else
			AzPML = 1.0;
	}
	else
	{
		if (k < KsMin)
			AzPML = pow(AttnH[KsMin - k - 1], -k0z);
		else if (k >= KsMax)
			AzPML = pow(AttnH[k - KsMax], k0z);
		else
			AzPML = 1.0;
	}

	APML = AxPML*AyPML*AzPML;
	return APML;
}