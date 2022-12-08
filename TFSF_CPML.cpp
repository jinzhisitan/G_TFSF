#include<cmath>
#include"TFSF.h"

void TFSF::add_TFSF_Box_E_CPML(const Matrix<double> &ce_x_1, const Matrix<double> &ce_x_2, const Matrix<double> &ce_y_1, 
	const Matrix<double> &ce_y_2, const Matrix<double> &ce_z_1, const Matrix<double> &ce_z_2,
	Matrix<double> &psi_Exy_1, Matrix<double> &psi_Exy_2, Matrix<double> &psi_Exz_1, Matrix<double> &psi_Exz_2, Matrix<double> &psi_Eyx_1, Matrix<double> &psi_Eyx_2, 
	Matrix<double> &psi_Eyz_1, Matrix<double> &psi_Eyz_2, Matrix<double> &psi_Ezx_1, Matrix<double> &psi_Ezx_2, Matrix<double> &psi_Ezy_1, Matrix<double> &psi_Ezy_2)
{
	add_TFSF_X1_E_CPML(ce_x_1, psi_Eyx_1, psi_Ezx_1);
	add_TFSF_X2_E_CPML(ce_x_2, psi_Eyx_2, psi_Ezx_2);
	add_TFSF_Y1_E_CPML(ce_y_1, psi_Exy_1, psi_Ezy_1);
	add_TFSF_Y2_E_CPML(ce_y_2, psi_Exy_2, psi_Ezy_2);
	add_TFSF_Z1_E_CPML(ce_z_1, psi_Exz_1, psi_Eyz_1);
	add_TFSF_Z2_E_CPML(ce_z_2, psi_Exz_2, psi_Eyz_2);
}

void TFSF::add_TFSF_Box_H_CPML(const Matrix<double> &ch_x_1, const Matrix<double> &ch_x_2, const Matrix<double> &ch_y_1,
	const Matrix<double> &ch_y_2, const Matrix<double> &ch_z_1, const Matrix<double> &ch_z_2,
	Matrix<double> &psi_Hxy_1, Matrix<double> &psi_Hxy_2, Matrix<double> &psi_Hxz_1, Matrix<double> &psi_Hxz_2, Matrix<double> &psi_Hyx_1, Matrix<double> &psi_Hyx_2,
	Matrix<double> &psi_Hyz_1, Matrix<double> &psi_Hyz_2, Matrix<double> &psi_Hzx_1, Matrix<double> &psi_Hzx_2, Matrix<double> &psi_Hzy_1, Matrix<double> &psi_Hzy_2)
{
	add_TFSF_X1_H_CPML(ch_x_1, psi_Hyx_1, psi_Hzx_1);
	add_TFSF_X2_H_CPML(ch_x_2, psi_Hyx_2, psi_Hzx_2);
	add_TFSF_Y1_H_CPML(ch_y_1, psi_Hxy_1, psi_Hzy_1);
	add_TFSF_Y2_H_CPML(ch_y_2, psi_Hxy_2, psi_Hzy_2);
	add_TFSF_Z1_H_CPML(ch_z_1, psi_Hxz_1, psi_Hyz_1);
	add_TFSF_Z2_H_CPML(ch_z_2, psi_Hxz_2, psi_Hyz_2);
}


////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////    E
void TFSF::add_TFSF_X1_E_CPML(const Matrix<double>&ce_x_1, Matrix<double> &psi_Eyx_1, Matrix<double> &psi_Ezx_1)
{
	////总场边界条件深入CPML
	if (ItMin > IsMin && ItMin<IsMax) {
		return;
	}

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

			APML = attenuationFactor(ItMin - 0.5, j + 0.5, k);
			psi_Eyx_1(ItMin, j, k) += ce_x_1(ItMin)*HzCB_x1*APML;   //表6-3 修改mu0

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
			psi_Ezx_1(ItMin, j, k) -= ce_x_1(ItMin) *HyCB_x1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_X2_E_CPML(const Matrix<double>&ce_x_2, Matrix<double> &psi_Eyx_2, Matrix<double> &psi_Ezx_2)
{
	////总场边界条件深入CPML
	if (ItMax < IsMax&&ItMax>IsMin) {
		return;
	}

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
			psi_Eyx_2(ItMax, j, k) -= ce_x_2(ItMax) *HzCB_x2*APML;   //表6-3 修改mu0
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
			psi_Ezx_2(ItMax, j, k) += ce_x_2(ItMax)  *HyCB_x2*APML;   //表6-3 修改mu0
		}
	}

}

void TFSF::add_TFSF_Y1_E_CPML(const Matrix<double>&ce_y_1, Matrix<double> &psi_Exy_1, Matrix<double> &psi_Ezy_1)
{
	////总场边界条件深入CPML
	if (JtMin > JsMin&&JtMin<JsMax) {
		return;
	}

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
			psi_Ezy_1(i, JtMin, k) += ce_y_1(JtMin)*HxCB_y1*APML;   //表6-3 修改mu0
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
			psi_Exy_1(i, JtMin, k) -=  ce_y_1(JtMin) *HzCB_y1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y2_E_CPML(const Matrix<double>&ce_y_2, Matrix<double> &psi_Exy_2, Matrix<double> &psi_Ezy_2)
{
	////总场边界条件深入CPML
	if (JtMax < JsMax&&JtMax>JsMin) {
		return;
	}

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
			psi_Ezy_2(i, JtMax, k) -= ce_y_2(JtMax)*HxCB_y2*APML;   //表6-3 修改mu0
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
			psi_Exy_2(i, JtMax, k) += ce_y_2(JtMax)  *HzCB_y2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z1_E_CPML(const Matrix<double>&ce_z_1, Matrix<double> &psi_Exz_1, Matrix<double> &psi_Eyz_1)
{
	////总场边界条件深入CPML
	if (KtMin > KsMin&&KtMin<KsMax) {
		return;
	}

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
			psi_Exz_1(i, j, KtMin) += ce_z_1(KtMin)*HyCB_z1*APML;   //表6-3 修改mu0
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

			APML = attenuationFactor(i, j + 0.5, KtMin - 0.5);
			psi_Eyz_1(i, j, KtMin) -= ce_z_1(KtMin) *HxCB_z1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z2_E_CPML(const Matrix<double>&ce_z_2, Matrix<double> &psi_Exz_2, Matrix<double> &psi_Eyz_2)
{
	////总场边界条件深入CPML
	if (KtMax < KsMax && KtMax>KsMin) {
		return;
	}

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
			psi_Exz_2(i, j, KtMax) -= ce_z_2(KtMax)*HyCB_z2*APML;   //表6-3 修改mu0
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
			psi_Eyz_2(i, j, KtMax) += ce_z_2(KtMax) *HxCB_z2*APML;   //表6-3 修改mu0
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////    H
void TFSF::add_TFSF_X1_H_CPML(const Matrix<double> &ch_x_1, Matrix<double> &psi_Hyx_1, Matrix<double> &psi_Hzx_1)
{
	////总场边界条件深入CPML
	if (ItMin > IsMin && ItMin <IsMax) {
		return;
	}

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
			psi_Hyx_1(ItMin - 1, j, k) -=  ch_x_1(ItMin - 1)*EzCB_x1*APML;   //表6-3 修改mu0
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

			APML = attenuationFactor(ItMin, j + 0.5, k);
			psi_Hzx_1(ItMin - 1, j, k) += ch_x_1(ItMin - 1) *EyCB_x1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_X2_H_CPML(const Matrix<double> &ch_x_2, Matrix<double> &psi_Hyx_2, Matrix<double> &psi_Hzx_2)
{
	////总场边界条件深入CPML
	if (ItMax < IsMax && ItMax > IsMin) {
		return;
	}

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
			psi_Hyx_2(ItMax, j, k) += ch_x_2(ItMax)*EzCB_x2*APML;   //表6-3 修改mu0
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
			psi_Hzx_2(ItMax, j, k) -=  ch_x_2(ItMax) *EyCB_x2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y1_H_CPML(const Matrix<double> &ch_y_1, Matrix<double> &psi_Hxy_1, Matrix<double> &psi_Hzy_1)
{
	////总场边界条件深入CPML
	if (JtMin > JsMin&&JtMin<JsMax) {
		return;
	}

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
			psi_Hzy_1(i, JtMin - 1, k) -= ch_y_1(JtMin - 1)*ExCB_y1*APML;   //表6-3 修改mu0
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
			psi_Hxy_1(i, JtMin - 1, k) += ch_y_1(JtMin - 1) *EzCB_y1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Y2_H_CPML(const Matrix<double> &ch_y_2, Matrix<double> &psi_Hxy_2, Matrix<double> &psi_Hzy_2)
{
	////总场边界条件深入CPML
	if (JtMax < JsMax && JtMax>JsMin) {
		return;
	}

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
			psi_Hzy_2(i, JtMax, k) +=ch_y_2(JtMax)*ExCB_y2*APML;   //表6-3 修改mu0
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
			psi_Hxy_2(i, JtMax, k) -=ch_y_2(JtMax) *EzCB_y2*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z1_H_CPML(const Matrix<double> &ch_z_1, Matrix<double> &psi_Hxz_1, Matrix<double> &psi_Hyz_1)
{
	////总场边界条件深入CPML
	if (KtMin > KsMin && KtMin<KsMax) {
		return;
	}

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
			psi_Hxz_1(i, j, KtMin - 1) -= ch_z_1(KtMin - 1)*EyCB_z1*APML;   //表6-3 修改mu0
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
			psi_Hyz_1(i, j, KtMin - 1) += ch_z_1(KtMin - 1) *ExCB_z1*APML;   //表6-3 修改mu0
		}
	}
}

void TFSF::add_TFSF_Z2_H_CPML(const Matrix<double> &ch_z_2, Matrix<double> &psi_Hxz_2, Matrix<double> &psi_Hyz_2)
{
	////总场边界条件深入CPML
	if (KtMax < KsMax && KtMax>KsMin) {
		return;
	}

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
			psi_Hxz_2(i, j, KtMax) +=  ch_z_2(KtMax)*EyCB_z2*APML;   //表6-3 修改mu0
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
			psi_Hyz_2(i, j, KtMax) -=  ch_z_2(KtMax)*ExCB_z2*APML;   //表6-3 修改mu0
		}
	}
}