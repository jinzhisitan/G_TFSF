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
	//Is = -40;
	return Is;
}


void TFSF::initializeTSFS()
{
	Ein.Allocate(IncidentStart, IncidentEnd);
	Hin.Allocate(IncidentStart, IncidentEnd - 1);
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


void TFSF::add_TFSF_Box_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ex, const Matrix<double> &den_ey,
	const Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez)
{
	add_TFSF_X1_E(CB, ob, den_ex, Ey, Ez);
	add_TFSF_X2_E(CB, ob, den_ex, Ey, Ez);
	add_TFSF_Y1_E(CB, ob, den_ey, Ex, Ez);
	add_TFSF_Y2_E(CB, ob, den_ey, Ex, Ez);
	add_TFSF_Z1_E(CB, ob, den_ez, Ex, Ey);
	add_TFSF_Z2_E(CB, ob, den_ez, Ex, Ey);
}

void TFSF::add_TFSF_X1_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ex, 
	Matrix<double> &Ey, Matrix<double> &Ez)
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
	double factor2;
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

			factor2 = CB(ob(ItMin, j, k),ob(ItMin - 1, j, k),ob(ItMin, j, k - 1),ob(ItMin - 1, j, k - 1));
			Ey(ItMin, j, k) += factor2*den_ex(ItMin)*HzCB_x1*APML;   //表6-3 修改mu0
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

			factor2 = CB(ob(ItMin, j, k),ob(ItMin - 1, j, k),ob(ItMin, j - 1, k),ob(ItMin - 1, j - 1, k));
			Ez(ItMin, j, k) -= factor2*den_ex(ItMin)*HyCB_x1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_X2_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ex, 
	Matrix<double> &Ey, Matrix<double> &Ez)
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
	double factor2;
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

			factor2 = CB(ob(ItMax, j, k), ob(ItMax - 1, j, k), ob(ItMax, j, k - 1), ob(ItMax - 1, j, k - 1));
			Ey(ItMax, j, k) -= factor2 *den_ex(ItMax)*HzCB_x2*APML;   //表6-3 修改mu0
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

			factor2 = CB(ob(ItMax, j, k),ob(ItMax - 1, j, k),ob(ItMax, j - 1, k),ob(ItMax - 1, j - 1, k));
			Ez(ItMax, j, k) += factor2*den_ex(ItMax) *HyCB_x2*APML;   //表6-3 修改mu0
		}
	}

}

void TFSF::add_TFSF_Y1_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ey,
	Matrix<double> &Ex, Matrix<double> &Ez)
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
	double factor2;
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

			factor2 = CB(ob(i, JtMin, k),ob(i - 1, JtMin, k),ob(i, JtMin - 1, k),ob(i - 1, JtMin - 1, k));
			Ez(i, JtMin, k) += factor2 *den_ey(JtMin)*HxCB_y1*APML;   //表6-3 修改mu0
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

			factor2 = CB(ob(i, JtMin, k),ob(i, JtMin, k - 1),ob(i, JtMin - 1, k),ob(i, JtMin - 1, k - 1));
			Ex(i, JtMin, k) -= factor2 *den_ey(JtMin)*HzCB_y1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y2_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ey, 
	Matrix<double> &Ex, Matrix<double> &Ez)
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
	double factor2;
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

			factor2 = CB(ob(i, JtMax, k),ob(i - 1, JtMax, k),ob(i, JtMax - 1, k),ob(i - 1, JtMax - 1, k));
			Ez(i, JtMax, k) -= factor2*den_ey(JtMax)*HxCB_y2*APML;   //表6-3 修改mu0
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

			factor2 = CB(ob(i, JtMax, k),ob(i, JtMax, k - 1),ob(i, JtMax - 1, k),ob(i, JtMax - 1, k - 1));
			Ex(i, JtMax, k) += factor2* den_ey(JtMax) *HzCB_y2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z1_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ez,
	Matrix<double> &Ex, Matrix<double> &Ey)
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
	double factor2;
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

			factor2 = CB(ob(i, j, KtMin), ob(i, j, KtMin - 1), ob(i, j - 1, KtMin), ob(i, j - 1, KtMin - 1));
			Ex(i, j, KtMin) += factor2*den_ez(KtMin)*HyCB_z1*APML;   //表6-3 修改mu0
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

			factor2 = CB(ob(i, j, KtMin),ob(i - 1, j, KtMin),ob(i, j, KtMin - 1),ob(i - 1, j, KtMin - 1));
			Ey(i, j, KtMin) -= factor2 *den_ez(KtMin)*HxCB_z1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z2_E(const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ez, 
	Matrix<double> &Ex, Matrix<double> &Ey)
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
	double factor2;
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

			factor2 = CB(ob(i, j, KtMax),ob(i, j, KtMax - 1),ob(i, j - 1, KtMax),ob(i, j - 1, KtMax - 1));
			Ex(i, j, KtMax) -= factor2*den_ez(KtMax)*HyCB_z2*APML;   //表6-3 修改mu0
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

			factor2 = CB(ob(i, j, KtMax),ob(i - 1, j, KtMax),ob(i, j, KtMax - 1),ob(i - 1, j, KtMax - 1));
			Ey(i, j, KtMax) += factor2*den_ez(KtMax)*HxCB_z2*APML;   //表6-3 修改mu0
		}
	}
}


void TFSF::add_TFSF_Box_H(const Matrix<double> &den_hx, const Matrix<double> &den_hy, const Matrix<double> &den_hz, 
	Matrix<double> &Hx, Matrix<double> &Hy, Matrix<double> &Hz)
{
	add_TFSF_X1_H(den_hx, Hy, Hz);
	add_TFSF_X2_H(den_hx, Hy, Hz);
	add_TFSF_Y1_H(den_hy, Hx, Hz);
	add_TFSF_Y2_H(den_hy, Hx, Hz);
	add_TFSF_Z1_H(den_hz, Hx, Hy);
	add_TFSF_Z2_H(den_hz, Hx, Hy);
}

void TFSF::add_TFSF_X1_H(const Matrix<double> &den_hx, Matrix<double> &Hy, Matrix<double> &Hz)
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
			Hy(ItMin - 1, j, k) -=  dt / mu0 *den_hx(ItMin - 1)*EzCB_x1*APML;   //表6-3 修改mu0
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
			Hz(ItMin - 1, j, k) += dt / mu0 *den_hx(ItMin - 1)*EyCB_x1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_X2_H(const Matrix<double> &den_hx, Matrix<double> &Hy, Matrix<double> &Hz)
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
			Hy(ItMax, j, k) += dt / mu0* den_hx(ItMax)*EzCB_x2*APML;   //表6-3 修改mu0
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
			Hz(ItMax, j, k) -= dt / mu0 * den_hx(ItMax)*EyCB_x2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y1_H(const Matrix<double> &den_hy, Matrix<double> &Hx, Matrix<double> &Hz)
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
			Hz(i, JtMin - 1, k) -= dt / mu0*den_hy(JtMin - 1)*ExCB_y1*APML;   //表6-3 修改mu0
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
			Hx(i, JtMin - 1, k) += dt / mu0*den_hy(JtMin - 1)*EzCB_y1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y2_H(const Matrix<double> &den_hy, Matrix<double> &Hx, Matrix<double> &Hz)
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
			Hz(i, JtMax, k) += dt / mu0*den_hy(JtMax)*ExCB_y2*APML;   //表6-3 修改mu0
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
			Hx(i, JtMax, k) -= dt / mu0*den_hy(JtMax)*EzCB_y2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z1_H(const Matrix<double> &den_hz, Matrix<double> &Hx, Matrix<double> &Hy)
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
			Hx(i, j, KtMin - 1) -= dt / mu0*den_hz(KtMin - 1)*EyCB_z1*APML;   //表6-3 修改mu0
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
			Hy(i, j, KtMin - 1) += dt / mu0*den_hz(KtMin - 1)*ExCB_z1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z2_H(const Matrix<double> &den_hz, Matrix<double> &Hx, Matrix<double> &Hy)
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
			Hx(i, j, KtMax) += dt / mu0* den_hz(KtMax)*EyCB_z2*APML;   //表6-3 修改mu0
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
			Hy(i, j, KtMax) -= dt / mu0* den_hz(KtMax)*ExCB_z2*APML;   //表6-3 修改mu0
		}
	}
}



void TFSF::add_TSFS_Box_E_analysis(double nt, const Matrix<double> &CB, const Matrix<int> &ob, const Matrix<double> &den_ex, const Matrix<double> &den_ey,
	const Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez)
{
	double APML = 1.0;
	double phase = 0.0;

	double kx = dx*sin(thi)*cos(phi);
	double ky = dy*sin(thi)*sin(phi);
	double kz = dz*cos(thi);
	double Hxtemp = -sin(phi)*cos(alpha) - cos(thi)*cos(phi)*sin(alpha);
	double Hytemp = cos(phi)*cos(alpha) - cos(thi)*sin(phi)*sin(alpha);
	double Hztemp = sin(thi)*sin(alpha);
	int i, j, k;
	double T1, tt;
	double factor2;

	//////////////////////////////////////////////////////////////////////////////////
	////   X面
	//////////////////////////////////////////////////////////////////////////////////
	double HzCB_x1, HzCB_x2;
	for (j = JtMin; j <= JtMax - 1; j++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = (ItMin - 0.5)*kx + (j + 0.5)*ky + k*kz;  //Hy
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(ItMin - 0.5, j + 0.5, k);
			phase = phaseCorrection(ItMin - 0.5, j + 0.5, k);
			factor2 = CB(ob(ItMin, j, k), ob(ItMin - 1, j, k), ob(ItMin, j, k - 1), ob(ItMin - 1, j, k - 1));
			HzCB_x1 = source(tt, phase) / Z0* Hztemp;
			Ey(ItMin, j, k) += factor2*den_ex(ItMin)*HzCB_x1*APML;   //表6-3 修改mu0

			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = (ItMax + 0.5)*kx + (j + 0.5)*ky + k*kz;  //Hy
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(ItMax + 0.5, j + 0.5, k);
			phase = phaseCorrection(ItMax + 0.5, j + 0.5, k);
			factor2 = CB(ob(ItMax, j, k), ob(ItMax - 1, j, k), ob(ItMax, j, k - 1), ob(ItMax - 1, j, k - 1));
			HzCB_x2 = source(tt, phase) / Z0* Hztemp;
			Ey(ItMax, j, k) -= factor2*den_ex(ItMax)*HzCB_x2*APML;   //表6-3 修改mu0
		}
	}
	//   x1  Ez
	double HyCB_x1, HyCB_x2;
	for (j = JtMin; j <= JtMax; j++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = (ItMin - 0.5)*kx + j *ky + (k + 0.5)*kz;   //Hx
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(ItMin - 0.5, j, k + 0.5);
			phase = phaseCorrection(ItMin - 0.5, j, k + 0.5);
			factor2 = CB(ob(ItMin, j, k), ob(ItMin - 1, j, k), ob(ItMin, j - 1, k), ob(ItMin - 1, j - 1, k));
			HyCB_x1 = source(tt, phase) / Z0*Hytemp;
			Ez(ItMin, j, k) -= factor2*den_ex(ItMin)*HyCB_x1*APML;   //表6-3 修改mu0

			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = (ItMax + 0.5)*kx + j *ky + (k + 0.5)*kz;   //Hx
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(ItMax + 0.5, j, k + 0.5);
			phase = phaseCorrection(ItMax + 0.5, j, k + 0.5);
			factor2 = CB(ob(ItMin, j, k), ob(ItMin - 1, j, k), ob(ItMin, j - 1, k), ob(ItMin - 1, j - 1, k));
			HyCB_x2 = source(tt, phase) / Z0*Hytemp;
			Ez(ItMax, j, k) += factor2*den_ex(ItMax)*HyCB_x2*APML;   //表6-3 修改mu0
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	////   Y面
	//////////////////////////////////////////////////////////////////////////////////
	double HxCB_y1, HxCB_y2;
	for (i = ItMin; i <= ItMax; i++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = i*kx + (JtMin - 0.5)*ky + (k + 0.5)*kz;  //Hx
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i, JtMin - 0.5, k + 0.5);
			phase = phaseCorrection(i, JtMin - 0.5, k + 0.5);
			factor2 = CB(ob(i, JtMin, k), ob(i - 1, JtMin, k), ob(i, JtMin - 1, k), ob(i - 1, JtMin - 1, k));
			HxCB_y1 = source(tt, phase) / Z0* Hxtemp;
			Ez(i, JtMin, k) += factor2 *den_ey(JtMin)*HxCB_y1*APML;   //表6-3 修改mu0


		    /////////////////////////////////////////////////////////////////////////////////////////
			T1 = i*kx + (JtMax + 0.5)*ky + (k + 0.5)*kz;  //Hx
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i, JtMax + 0.5, k + 0.5);
			phase = phaseCorrection(i, JtMax + 0.5, k + 0.5);
			factor2 = CB(ob(i, JtMax, k), ob(i - 1, JtMax, k), ob(i, JtMax - 1, k), ob(i - 1, JtMax - 1, k));
			HxCB_y2 = source(tt, phase) / Z0* Hxtemp;
			Ez(i, JtMax, k) -= factor2 *den_ey(JtMax)*HxCB_y2*APML;   //表6-3 修改mu0
		}
	}
	//    y1  Ex
	double HzCB_y1, HzCB_y2;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = (i + 0.5)*kx + (JtMin - 0.5)*ky + k *kz;  //Hz
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i + 0.5, JtMin - 0.5, k);
			phase = phaseCorrection(i + 0.5, JtMin - 0.5, k);
			factor2 = CB(ob(i, JtMin, k), ob(i, JtMin, k - 1), ob(i, JtMin - 1, k), ob(i, JtMin - 1, k - 1));
			HzCB_y1 = source(tt, phase) / Z0*Hztemp;
			Ex(i, JtMin, k) -= factor2 *den_ey(JtMin)*HzCB_y1*APML;   //表6-3 修改mu0

			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = (i + 0.5)*kx + (JtMax + 0.5)*ky + k *kz;  //Hz
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i + 0.5, JtMax + 0.5, k);
			phase = phaseCorrection(i + 0.5, JtMax + 0.5, k);
			factor2 = CB(ob(i, JtMax, k), ob(i, JtMax, k - 1), ob(i, JtMax - 1, k), ob(i, JtMax - 1, k - 1));
			HzCB_y2 = source(tt, phase) / Z0*Hztemp;
			Ex(i, JtMax, k) += factor2* den_ey(JtMax)*HzCB_y2*APML;   //表6-3 修改mu0
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	////   Z面
	//////////////////////////////////////////////////////////////////////////////////
	double HyCB_z1, HyCB_z2;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (j = JtMin; j <= JtMax; j++) {
			T1 = (i + 0.5)*kx + j*ky + (KtMin - 0.5)*kz;  //Hy
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i + 0.5, j, KtMin - 0.5);
			phase = phaseCorrection(i + 0.5, j, KtMin - 0.5);
			factor2 = CB(ob(i, j, KtMin), ob(i, j, KtMin - 1), ob(i, j - 1, KtMin), ob(i, j - 1, KtMin - 1));
			HyCB_z1 = source(tt, phase) / Z0* Hytemp;
			Ex(i, j, KtMin) += factor2*den_ez(KtMin)*HyCB_z1*APML;   //表6-3 修改mu0


			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = (i + 0.5)*kx + j*ky + (KtMax + 0.5)*kz;  //Hy
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i + 0.5, j, KtMax + 0.5);
			phase = phaseCorrection(i + 0.5, j, KtMax + 0.5);
			factor2 = CB(ob(i, j, KtMax), ob(i, j, KtMax - 1), ob(i, j - 1, KtMax), ob(i, j - 1, KtMax - 1));
			HyCB_z2 = source(tt, phase) / Z0*Hytemp;
			Ex(i, j, KtMax) -= factor2*den_ez(KtMax)*HyCB_z2*APML;   //表6-3 修改mu0
		}
	}
	//   z2  Ey
	double HxCB_z1, HxCB_z2;
	for (i = ItMin; i <= ItMax; i++) {
		for (j = JtMin; j <= JtMax - 1; j++) {
			T1 = i*kx + (j + 0.5)*ky + (KtMin - 0.5)*kz;   //Hx
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i, j + 0.5, KtMin - 0.5);
			phase = phaseCorrection(i, j + 0.5, KtMin - 0.5);
			factor2 = CB(ob(i, j, KtMin), ob(i - 1, j, KtMin), ob(i, j, KtMin - 1), ob(i - 1, j, KtMin - 1));
			HxCB_z1 = source(tt, phase) / Z0*Hxtemp;
			Ey(i, j, KtMin) -= factor2 *den_ez(KtMin)*HxCB_z1*APML;   //表6-3 修改mu0


			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = i*kx + (j + 0.5)*ky + (KtMax + 0.5)*kz;   //Hx
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i, j + 0.5, KtMax + 0.5);
			phase = phaseCorrection(i, j + 0.5, KtMax + 0.5);
			factor2 = CB(ob(i, j, KtMax), ob(i - 1, j, KtMax), ob(i, j, KtMax - 1), ob(i - 1, j, KtMax - 1));
			HxCB_z2 = source(tt, phase) / Z0*Hxtemp;
			Ey(i, j, KtMax) += factor2*den_ez(KtMax)*HxCB_z2*APML;   //表6-3 修改mu0
		}
	}
}



void TFSF::add_TSFS_Box_H_analysis(double nt, const Matrix<double> &den_hx, const Matrix<double> &den_hy, const Matrix<double> &den_hz, 
	Matrix<double> &Hx, Matrix<double> &Hy, Matrix<double> &Hz)
{
	double APML = 1.0;
	double phase = 0.0;


	double kx = dx*sin(thi)*cos(phi);
	double ky = dy*sin(thi)*sin(phi);
	double kz = dz*cos(thi);
	double Extemp = -sin(phi)*sin(alpha) + cos(thi)*cos(phi)*cos(alpha);
	double Eytemp = cos(phi)*sin(alpha) + cos(thi)*sin(phi)*cos(alpha);
	double Eztemp = -sin(thi)*cos(alpha);
	int i, j, k;
	double T1, tt;

	//////////////////////////////////////////////////////////////////////////////////
	////   X面
	//////////////////////////////////////////////////////////////////////////////////
	//x1 Hy
	double EzCB_x1, EzCB_x2;
	for (j = JtMin; j <= JtMax; j++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = ItMin*kx + j*ky + (k + 0.5)*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(ItMin, j, k + 0.5);
			phase = phaseCorrection(ItMin, j, k + 0.5);
			EzCB_x1 = source(tt, phase)* Eztemp;
			Hy(ItMin - 1, j, k) -= dt / mu0 *den_hx(ItMin - 1)*EzCB_x1*APML;   //表6-3 修改mu0


			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = ItMax*kx + j*ky + (k + 0.5)*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(ItMax, j, k + 0.5);
			phase = phaseCorrection(ItMax, j, k + 0.5);
			EzCB_x2 = source(tt, phase)* Eztemp;
			Hy(ItMax, j, k) += dt / mu0 * den_hx(ItMax)*EzCB_x2*APML;   //表6-3 修改mu0
		}
	}
	//x1 Hz
	double EyCB_x1, EyCB_x2;
	for (j = JtMin; j <= JtMax - 1; j++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = ItMin*kx + (j + 0.5)*ky + k *kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(ItMin, j + 0.5, k);
			phase = phaseCorrection(ItMin, j + 0.5, k);
			EyCB_x1 = source(tt, phase)*Eytemp;
			Hz(ItMin - 1, j, k) += dt / mu0 *den_hx(ItMin - 1)*EyCB_x1*APML;   //表6-3 修改mu0


			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = ItMax*kx + (j + 0.5)*ky + k *kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(ItMax, j + 0.5, k);
			phase = phaseCorrection(ItMax, j + 0.5, k);
			EyCB_x2 = source(tt, phase)*Eytemp;
			Hz(ItMax, j, k) -= dt / mu0 *den_hx(ItMax)*EyCB_x2*APML;   //表6-3 修改mu0
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	////   Y面
	//////////////////////////////////////////////////////////////////////////////////
	double ExCB_y1, ExCB_y2;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (k = KtMin; k <= KtMax; k++) {
			T1 = (i + 0.5)*kx + JtMin*ky + k*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i + 0.5, JtMin, k);
			phase = phaseCorrection(i + 0.5, JtMin, k);
			ExCB_y1 = source(tt, phase)* Extemp;
			Hz(i, JtMin - 1, k) -= dt / mu0 *den_hy(JtMin - 1)*ExCB_y1*APML;   //表6-3 修改mu0

			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = (i + 0.5)*kx + JtMax*ky + k*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i + 0.5, JtMax, k);
			phase = phaseCorrection(i + 0.5, JtMax, k);
			ExCB_y2 = source(tt, phase)* Extemp;
			Hz(i, JtMax, k) += dt / mu0 *den_hy(JtMax)*ExCB_y2*APML;   //表6-3 修改mu0
		}
	}
	// y1  Hx
	double EzCB_y1, EzCB_y2;
	for (i = ItMin; i <= ItMax; i++) {
		for (k = KtMin; k <= KtMax - 1; k++) {
			T1 = i*kx + JtMin*ky + (k + 0.5)*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i, JtMin, k + 0.5);
			phase = phaseCorrection(i, JtMin, k + 0.5);
			EzCB_y1 = source(tt, phase)* Eztemp;
			Hx(i, JtMin - 1, k) += dt / mu0 *den_hy(JtMin - 1)*EzCB_y1*APML;   //表6-3 修改mu0


			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = i*kx + JtMax*ky + (k + 0.5)*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i, JtMax, k + 0.5);
			phase = phaseCorrection(i, JtMax, k + 0.5);
			EzCB_y2 = source(tt, phase)* Eztemp;
			Hx(i, JtMax, k) -= dt / mu0 *den_hy(JtMax)*EzCB_y2*APML;   //表6-3 修改mu0
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	////   Z面
	//////////////////////////////////////////////////////////////////////////////////
	double EyCB_z1, EyCB_z2;
	for (i = ItMin; i <= ItMax; i++) {
		for (j = JtMin; j <= JtMax - 1; j++) {
			T1 = i*kx + (j + 0.5)*ky + KtMin*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i, j + 0.5, KtMin);
			phase = phaseCorrection(i, j + 0.5, KtMin);
			EyCB_z1 = source(tt, phase)*Eytemp;
			Hx(i, j, KtMin - 1) -= dt / mu0 *den_hz(KtMin - 1)*EyCB_z1*APML;   //表6-3 修改mu0


			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = i*kx + (j + 0.5)*ky + KtMax*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i, j + 0.5, KtMax);
			phase = phaseCorrection(i, j + 0.5, KtMax);
			EyCB_z2 = source(tt, phase)*Eytemp;
			Hx(i, j, KtMax) += dt / mu0 * den_hz(KtMax)*EyCB_z2*APML;   //表6-3 修改mu0
		}
	}

	double ExCB_z1, ExCB_z2;
	for (i = ItMin; i <= ItMax - 1; i++) {
		for (j = JtMin; j <= JtMax; j++) {
			T1 = (i + 0.5)*kx + j*ky + KtMin*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i + 0.5, j, KtMin);
			phase = phaseCorrection(i + 0.5, j, KtMin);
			ExCB_z1 = source(tt, phase)*Extemp;
			Hy(i, j, KtMin - 1) += dt / mu0 *den_hz(KtMin - 1)*ExCB_z1*APML;   //表6-3 修改mu0
			
			/////////////////////////////////////////////////////////////////////////////////////////
			T1 = (i + 0.5)*kx + j*ky + KtMax*kz;
			T1 -= Isource*delta;
			tt = nt*dt - T1 / c0;
			APML = attenuationFactor(i + 0.5, j, KtMax);
			phase = phaseCorrection(i + 0.5, j, KtMax);
			ExCB_z2 = source(tt, phase)* Extemp;
			Hy(i, j, KtMax) -= dt / mu0 * den_hz(KtMax)*ExCB_z2*APML;   //表6-3 修改mu0
		}
	}
}