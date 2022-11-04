#include<iomanip>
#include"TFSF.h"

std::string fname = "zz_TFSF.dat";
ofstream fout(fname);

TFSF::TFSF()
{
	std::cout << "TFSF: void constructor\n";
}


TFSF::TFSF(double dx_, double dy_, double dz_, double alpha_, double thi_, double phi_, double dt_, int IncidentStart_, int IncidentEnd_,int ItMin_, 
	int ItMax_, int JtMin_, int JtMax_, int KtMin_, int KtMax_, int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_,int KsMax_)
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

	Ein.Allocate(IncidentStart, IncidentEnd);
	Hin.Allocate(IncidentStart, IncidentEnd - 1);
	APML();
}

double TFSF::ComputeComponentFrac()
{
	//PURPOSE:
	//	计算一维辅助入射波空间离散间隔与三维离散间隔的比值：Delta/dx
	//Reference:
	//	[1]  葛德彪, 闫玉波. (2011). 电磁波时域有限差分方法（第三版）.西安电子科技大学出版社.
	//  [2]  Taflove, A., &Hagness, S.C. (2005).Computational electrodynamics 
	//       : the finite - difference time - domain method 3rd ed. Artech house.
	double st = sin(thi);
	double ct = cos(thi);
	double sp = sin(phi);
	double cp = cos(phi);

	//ref[2], P217, Eq(5.69)
	double m_VFrac = sqrt(pow(st, 4)*(pow(cp, 4) + pow(sp, 4)) + pow(ct, 4));

	//m_VFrac = 1.0;
	std::cout << "m_VFrac= " << m_VFrac << endl;
	return m_VFrac;
}

int TFSF::source_position()
{
	const double pi = 3.14159265358979323846;
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





void TFSF::update1D_Einc(double nt)
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
	Ein(IncidentStart) = EBin[1] + (c0*dt - delta) / (c0*dt + delta)*(Ein(IncidentStart + 1) - EBin[0]);    //一维入射波吸收边界 Eq(4-3-17)
	Ein(IncidentEnd) = EBin[2] + (c0*dt - delta) / (c0*dt + delta)*(Ein(IncidentEnd - 1) - EBin[3]);
	//Ein(IncidentStart) = EBin[0] - c0*dt / delta*(EBin[0] - EBin[1]);    //一维入射波吸收边界 Eq(4-3-20)
	//Ein(IncidentEnd) = EBin[3] - c0*dt / delta*(EBin[3] - EBin[2]);
	Ein(Isource) = source(nt*dt);

	fout << scientific;
	fout.precision(4);
	fout << setw(11) << nt*dt*1E6 << setw(14) << Ein(Isource)
		<< setw(14) << Ein(0) << endl;
}

void TFSF::update1D_Hinc()
{
	double FH = dt / mu0 / delta;
	// 产生入射波  Ein  Hin
	for (int i = IncidentStart; i <= IncidentEnd - 1; i++) {
		Hin(i) = Hin(i) - FH*(Ein(i + 1) - Ein(i));
	}
}

void TFSF::add_TFSF_Box_H(Matrix<double> &den_hx, Matrix<double> &den_hy, Matrix<double> &den_hz, 
	Matrix<double> &Hx, Matrix<double> &Hy, Matrix<double> &Hz)
{
	add_TFSF_X1_H(den_hx, Hy, Hz);
	add_TFSF_X2_H(den_hx, Hy, Hz);
	add_TFSF_Y1_H(den_hy, Hx, Hz);
	add_TFSF_Y2_H(den_hy, Hx, Hz);
	add_TFSF_Z1_H(den_hz, Hx, Hy);
	add_TFSF_Z2_H(den_hz, Hx, Hy);
}

void TFSF::add_TFSF_X1_H(Matrix<double> &den_hx, Matrix<double> &Hy, Matrix<double> &Hz)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	if (ItMin <= IsMin)
		AxPML = pow(AttnE[IsMin - ItMin], -k0x);
	else
		AxPML = 1.0;

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

			//******** y 方向
			if (j < JsMin)
				AyPML = pow(AttnE[JsMin - j], -k0y);
			else if (j > JsMax)
				AyPML = pow(AttnE[j - JsMax], k0y);
			else
				AyPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnH[KsMin - k - 1], -k0z);
			else if (k >= KsMax)
				AzPML = pow(AttnH[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

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

			//******** y 方向
			if (j < JsMin)
				AyPML = pow(AttnH[JsMin - j - 1], -k0y);
			else if (j >= JsMax)
				AyPML = pow(AttnH[j - JsMax], k0y);
			else
				AyPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnE[KsMin - k], -k0z);
			else if (k > KsMax)
				AzPML = pow(AttnE[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			Hz(ItMin - 1, j, k) += dt / mu0 *den_hx(ItMin - 1)*EyCB_x1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_X2_H(Matrix<double> &den_hx, Matrix<double> &Hy, Matrix<double> &Hz)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	if (ItMax >= IsMax)
		AxPML = pow(AttnE[ItMax - IsMax], k0x);
	else
		AxPML = 1.0;

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

			//******** y 方向
			if (j < JsMin)
				AyPML = pow(AttnE[JsMin - j], -k0y);
			else if (j > JsMax)
				AyPML = pow(AttnE[j - JsMax], k0y);
			else
				AyPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnH[KsMin - k - 1], -k0z);
			else if (k >= KsMax)
				AzPML = pow(AttnH[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			Hy(ItMax, j, k) +=  dt / mu0* den_hx(ItMax)*EzCB_x2*APML;   //表6-3 修改mu0
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

			//******** y 方向
			if (j < JsMin)
				AyPML = pow(AttnH[JsMin - j - 1], -k0y);
			else if (j >= JsMax)
				AyPML = pow(AttnH[j - JsMax], k0y);
			else
				AyPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnE[KsMin - k], -k0z);
			else if (k > KsMax)
				AzPML = pow(AttnE[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			Hz(ItMax, j, k) -=  dt / mu0 * den_hx(ItMax)*EyCB_x2;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y1_H(Matrix<double> &den_hy, Matrix<double> &Hx, Matrix<double> &Hz)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	if (JtMin <= JsMin)
		AyPML = pow(AttnE[JsMin - JtMin], -k0y);
	else
		AyPML = 1.0;

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

			//******** x 方向
			if (i < IsMin)
				AxPML = pow(AttnH[IsMin - i - 1], -k0x);
			else if (i >= IsMax)
				AxPML = pow(AttnH[i - IsMax], k0x);
			else
				AxPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnE[KsMin - k], -k0z);
			else if (k > KsMax)
				AzPML = pow(AttnE[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			Hz(i, JtMin - 1, k) -=  dt / mu0*den_hy(JtMin - 1)*ExCB_y1*APML;   //表6-3 修改mu0
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

			//******** y 方向
			if (i < IsMin)
				AxPML = pow(AttnE[IsMin - i], -k0x);
			else if (i > IsMax)
				AxPML = pow(AttnE[i - IsMax], k0x);
			else
				AxPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnH[KsMin - k - 1], -k0z);
			else if (k >= KsMax)
				AzPML = pow(AttnH[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			Hx(i, JtMin - 1, k) +=  dt / mu0*den_hy(JtMin - 1)*EzCB_y1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y2_H(Matrix<double> &den_hy, Matrix<double> &Hx, Matrix<double> &Hz)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	if (JtMax >= JsMax)
		AyPML = pow(AttnE[JtMax - JsMax], k0y);
	else
		AyPML = 1.0;

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

			//******** x 方向
			if (i < IsMin)
				AxPML = pow(AttnH[IsMin - i - 1], -k0x);
			else if (i >= IsMax)
				AxPML = pow(AttnH[i - IsMax], k0x);
			else
				AxPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnE[KsMin - k], -k0z);
			else if (k > KsMax)
				AzPML = pow(AttnE[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

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

			//******** x 方向
			if (i < IsMin)
				AxPML = pow(AttnE[IsMin - i], -k0x);
			else if (i > IsMax)
				AxPML = pow(AttnE[i - IsMax], k0x);
			else
				AxPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnH[KsMin - k - 1], -k0z);
			else if (k >= KsMax)
				AzPML = pow(AttnH[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			Hx(i, JtMax, k) -=  dt / mu0*den_hy(JtMax)*EzCB_y2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z1_H(Matrix<double> &den_hz, Matrix<double> &Hx, Matrix<double> &Hy)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	if (KtMin <= KsMin)
		AzPML = pow(AttnE[KsMin - KtMin], -k0z);
	else
		AzPML = 1.0;

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

			/////////////////////////////////////////////
			if (i < IsMin)
				AxPML = pow(AttnE[IsMin - i], -k0x);
			else if (i > IsMax)
				AxPML = pow(AttnE[i - IsMax], k0x);
			else
				AxPML = 1.0;

			if (j < JsMin)
				AyPML = pow(AttnH[JsMin - j - 1], -k0y);
			else if (j >= JsMax)
				AyPML = pow(AttnH[j - JsMax], k0y);
			else
				AyPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

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

			////////////////////////////////////////////////
			if (i < IsMin)
				AxPML = pow(AttnH[IsMin - i - 1], -k0x);
			else if (i >= IsMax)
				AxPML = pow(AttnH[i - IsMax], k0x);
			else
				AxPML = 1.0;

			if (j < JsMin)
				AyPML = pow(AttnE[JsMin - j], -k0y);
			else if (j > JsMax)
				AyPML = pow(AttnE[j - JsMax], k0y);
			else
				AyPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			Hy(i, j, KtMin - 1) += dt / mu0*den_hz(KtMin - 1)*ExCB_z1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z2_H(Matrix<double> &den_hz, Matrix<double> &Hx, Matrix<double> &Hy)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	
	if (KtMax >= KsMax)
		AzPML = pow(AttnE[KtMax - KsMax], k0z);
	else
		AzPML = 1.0;

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

			/////////////////////////////////////////////
			if (i < IsMin)
				AxPML = pow(AttnE[IsMin - i], -k0x);
			else if (i > IsMax)
				AxPML = pow(AttnE[i - IsMax], k0x);
			else
				AxPML = 1.0;

			if (j < JsMin)
				AyPML = pow(AttnH[JsMin - j - 1], -k0y);
			else if (j >= JsMax)
				AyPML = pow(AttnH[j - JsMax], k0y);
			else
				AyPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

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

			////////////////////////////////////////////////
			if (i < IsMin)
				AxPML = pow(AttnH[IsMin - i - 1], -k0x);
			else if (i >= IsMax)
				AxPML = pow(AttnH[i - IsMax], k0x);
			else
				AxPML = 1.0;

			if (j < JsMin)
				AyPML = pow(AttnE[JsMin - j], -k0y);
			else if (j > JsMax)
				AyPML = pow(AttnE[j - JsMax], k0y);
			else
				AyPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			Hy(i, j, KtMax) -= dt / mu0* den_hz(KtMax)*ExCB_z2*APML;   //表6-3 修改mu0
		}
	}
}


void TFSF::add_TFSF_Box_E(double CB[][MediaNo][MediaNo][MediaNo],Matrix<int> &ob, 
	Matrix<double> &den_ex, Matrix<double> &den_ey, Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey, Matrix<double> &Ez)
{
	add_TFSF_X1_E(CB, ob, den_ex, Ey, Ez);
	add_TFSF_X2_E(CB, ob, den_ex, Ey, Ez);
	add_TFSF_Y1_E(CB, ob, den_ey, Ex, Ez);
	add_TFSF_Y2_E(CB, ob, den_ey, Ex, Ez);
	add_TFSF_Z1_E(CB, ob, den_ez, Ex, Ey);
	add_TFSF_Z2_E(CB, ob, den_ez, Ex, Ey);
}

void TFSF::add_TFSF_X1_E(double CB[][MediaNo][MediaNo][MediaNo], 
	Matrix<int> &ob, Matrix<double> &den_ex, Matrix<double> &Ey, Matrix<double> &Ez)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	if (ItMin <= IsMin)
		AxPML = pow(AttnH[IsMin - ItMin], -k0x);
	else
		AxPML = 1.0;


	//入射角phi=180  thi=90-180
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

			//******** y 方向
			if (j < JsMin) 
				AyPML = pow(AttnH[JsMin - j - 1], -k0y);
			else if (j>=JsMax)
				AyPML = pow(AttnH[j - JsMax], k0y);
			else
				AyPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnE[KsMin - k], -k0z);
			else if (k > KsMax)
				AzPML = pow(AttnE[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(ItMin, j, k)][ob(ItMin - 1, j, k)][ob(ItMin, j, k - 1)][ob(ItMin - 1, j, k - 1)];
			Ey(ItMin, j, k) +=  factor2*den_ex(ItMin)*HzCB_x1*APML;   //表6-3 修改mu0
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

			//******** y 方向
			if (j < JsMin)
				AyPML = pow(AttnE[JsMin - j], -k0y);
			else if (j > JsMax)
				AyPML = pow(AttnE[j - JsMax], k0y);
			else
				AyPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnH[KsMin - k-1], -k0z);
			else if (k >= KsMax)
				AzPML = pow(AttnH[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(ItMin, j, k)][ob(ItMin-1, j, k)][ob(ItMin, j - 1, k)][ob(ItMin-1, j - 1, k)];
			Ez(ItMin, j, k) -= factor2*den_ex(ItMin)*HyCB_x1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_X2_E(double CB[][MediaNo][MediaNo][MediaNo],
	Matrix<int> &ob, Matrix<double> &den_ex, Matrix<double> &Ey, Matrix<double> &Ez)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML,AyPML, AzPML, APML;
	if (ItMax >= IsMax)
		AxPML = pow(AttnH[ItMax - IsMax], k0x);
	else
		AxPML = 1.0;

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

			//******** y 方向
			if (j < JsMin)
				AyPML = pow(AttnH[JsMin - j - 1], -k0y);
			else if (j >= JsMax)
				AyPML = pow(AttnH[j - JsMax], k0y);
			else
				AyPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnE[KsMin - k], -k0z);
			else if (k > KsMax)
				AzPML = pow(AttnE[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(ItMax, j, k)][ob(ItMax - 1, j, k)][ob(ItMax, j, k - 1)][ob(ItMax - 1, j, k - 1)];
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

			//******** y 方向
			if (j < JsMin)
				AyPML = pow(AttnE[JsMin - j], -k0y);
			else if (j > JsMax)
				AyPML = pow(AttnE[j - JsMax], k0y);
			else
				AyPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnH[KsMin - k - 1], -k0z);
			else if (k >= KsMax)
				AzPML = pow(AttnH[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(ItMax, j, k)][ob(ItMax-1, j, k)][ob(ItMax, j - 1, k)][ob(ItMax-1, j - 1, k)];
			Ez(ItMax, j, k) += factor2*den_ex(ItMax) *HyCB_x2*APML;   //表6-3 修改mu0
		}
	}

}

void TFSF::add_TFSF_Y1_E(double CB[][MediaNo][MediaNo][MediaNo],
	Matrix<int> &ob, Matrix<double> &den_ey, Matrix<double> &Ex, Matrix<double> &Ez)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	if (JtMin <= JsMin)
		AyPML = pow(AttnH[JsMin - JtMin], -k0y);
	else
		AyPML = 1.0;

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

			//******** x 方向
			if (i < IsMin)
				AxPML = pow(AttnE[IsMin - i], -k0x);
			else if (i > IsMax)
				AxPML = pow(AttnE[i - IsMax], k0x);
			else
				AxPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnH[KsMin - k-1], -k0z);
			else if (k >= KsMax)
				AzPML = pow(AttnH[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(i, JtMin, k)][ob(i - 1, JtMin, k)][ob(i, JtMin - 1, k)][ob(i - 1, JtMin - 1, k)];
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

			//******** x 方向
			if (i < IsMin)
				AxPML = pow(AttnH[IsMin - i-1], -k0x);
			else if (i >= IsMax)
				AxPML = pow(AttnH[i - IsMax], k0x);
			else
				AxPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnE[KsMin - k], -k0z);
			else if (k > KsMax)
				AzPML = pow(AttnE[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(i, JtMin, k)][ob(i, JtMin, k - 1)][ob(i, JtMin - 1, k)][ob(i, JtMin - 1, k - 1)];
			Ex(i, JtMin, k) -= factor2 *den_ey(JtMin)*HzCB_y1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y2_E(double CB[][MediaNo][MediaNo][MediaNo],
	Matrix<int> &ob, Matrix<double> &den_ey, Matrix<double> &Ex, Matrix<double> &Ez)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	if (JtMax >= JsMax)
		AyPML = pow(AttnH[JtMax - JsMax], k0y);
	else
		AyPML = 1.0;


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

			//******** x 方向
			if (i < IsMin)
				AxPML = pow(AttnE[IsMin - i], -k0x);
			else if (i > IsMax)
				AxPML = pow(AttnE[i - IsMax], k0x);
			else
				AxPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnH[KsMin - k - 1], -k0z);
			else if (k >= KsMax)
				AzPML = pow(AttnH[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(i, JtMax, k)][ob(i - 1, JtMax, k)][ob(i, JtMax - 1, k)][ob(i - 1, JtMax - 1, k)];
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

			//******** x 方向
			if (i < IsMin)
				AxPML = pow(AttnH[IsMin - i - 1], -k0x);
			else if (i >= IsMax)
				AxPML = pow(AttnH[i - IsMax], k0x);
			else
				AxPML = 1.0;
			//******** z 方向
			if (k < KsMin)
				AzPML = pow(AttnE[KsMin - k], -k0z);
			else if (k > KsMax)
				AzPML = pow(AttnE[k - KsMax], k0z);
			else
				AzPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(i, JtMax, k)][ob(i, JtMax, k - 1)][ob(i, JtMax - 1, k)][ob(i, JtMax - 1, k - 1)];
			Ex(i, JtMax, k) += factor2* den_ey(JtMax) *HzCB_y2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z1_E(double CB[][MediaNo][MediaNo][MediaNo],
	Matrix<int> &ob, Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;
	if (KtMin <= KsMin)
		AzPML = pow(AttnH[KsMin - KtMin], -k0z);
	else
		AzPML = 1.0;


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
			
			/////////////////////////////////////
			if (i < IsMin)
				AxPML = pow(AttnH[IsMin - i - 1], -k0x);
			else if (i >= IsMax)
				AxPML = pow(AttnH[i - IsMax], k0x);
			else
				AxPML = 1.0;

			if (j < JsMin)
				AyPML = pow(AttnE[JsMin - j], -k0y);
			else if (j > JsMax)
				AyPML = pow(AttnE[j - JsMax], k0y);
			else
				AyPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(i, j, KtMin)][ob(i, j, KtMin - 1)][ob(i, j - 1, KtMin)][ob(i, j - 1, KtMin - 1)];
			Ex(i, j, KtMin) += factor2*den_ez(KtMin)*HyCB_z1*APML;   //表6-3 修改mu0
		}
	}
	//   z1  Ey
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

			if (i < IsMin)
				AxPML = pow(AttnE[IsMin - i], -k0x);
			else if (i > IsMax)
				AxPML = pow(AttnE[i - IsMax], k0x);
			else
				AxPML = 1.0;

			if (j < JsMin)
				AyPML = pow(AttnH[JsMin - j-1], -k0y);
			else if (j >= JsMax)
				AyPML = pow(AttnH[j - JsMax], k0y);
			else
				AyPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(i, j, KtMin)][ob(i - 1, j, KtMin)][ob(i, j, KtMin - 1)][ob(i - 1, j, KtMin - 1)];
			Ey(i, j, KtMin) -= factor2 *den_ez(KtMin)*HxCB_z1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z2_E(double CB[][MediaNo][MediaNo][MediaNo],
	Matrix<int> &ob, Matrix<double> &den_ez, Matrix<double> &Ex, Matrix<double> &Ey)
{
	////总场边界条件深入CPML
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double AxPML, AyPML, AzPML, APML;

	if (KtMax >= KsMax)
		AzPML = pow(AttnH[KtMax - KsMax], k0z);
	else
		AzPML = 1.0;

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

			////////////////////////////////////////////
			if (i < IsMin)
				AxPML = pow(AttnH[IsMin - i - 1], -k0x);
			else if (i >= IsMax)
				AxPML = pow(AttnH[i - IsMax], k0x);
			else
				AxPML = 1.0;

			if (j < JsMin)
				AyPML = pow(AttnE[JsMin - j], -k0y);
			else if (j > JsMax)
				AyPML = pow(AttnE[j - JsMax], k0y);
			else
				AyPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(i, j, KtMax)][ob(i, j, KtMax - 1)][ob(i, j - 1, KtMax)][ob(i, j - 1, KtMax - 1)];
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

			/////////////////////////////////////////////////
			if (i < IsMin)
				AxPML = pow(AttnE[IsMin - i], -k0x);
			else if (i > IsMax)
				AxPML = pow(AttnE[i - IsMax], k0x);
			else
				AxPML = 1.0;

			if (j < JsMin)
				AyPML = pow(AttnH[JsMin - j - 1], -k0y);
			else if (j >= JsMax)
				AyPML = pow(AttnH[j - JsMax], k0y);
			else
				AyPML = 1.0;
			////////////////////////////////////////////////
			APML = AxPML*AyPML*AzPML;

			factor2 = CB[ob(i, j, KtMax)][ob(i - 1, j, KtMax)][ob(i, j, KtMax - 1)][ob(i - 1, j, KtMax - 1)];
			Ey(i, j, KtMax) += factor2*den_ez(KtMax)*HxCB_z2*APML;   //表6-3 修改mu0
		}
	}
}

