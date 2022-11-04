#include <string>
#include <fstream>
#include <iomanip>
#include "FDTD.h"

void FDTD::compute()
{
	double tempx, tempy, tempz, temp;
	Matrix<double> Emx(Imin, Imax, Jmin, Jmax, Kmin, Kmax);
	Matrix<double> Emy(Imin, Imax, Jmin, Jmax, Kmin, Kmax);
	Matrix<double> Emz(Imin, Imax, Jmin, Jmax, Kmin, Kmax);

	string f_point_Ex = "zz_point_Ex.dat";
	string f_point_Ey = "zz_point_Ey.dat";
	string f_point_Ez = "zz_point_Ez.dat";

	string f_line_E = "zz_line_E.dat";

	string f_Yface_Ex = "zz_Y=cons_Ex.dat";
	string f_Yface_Ey = "zz_Y=cons_Ey.dat";
	string f_Yface_Ez = "zz_Y=cons_Ez.dat";
	//string f_Yface_Et = "zz_Y=cons_Et.dat";

	string f_Zface_Ex = "zz_Z=cons_Ex.dat";
	string f_Zface_Ey = "zz_Z=cons_Ey.dat";
	string f_Zface_Ez = "zz_Z=cons_Ez.dat";
	//string f_Zface_Et = "zz_Z=cons_Et.dat";

	ofstream fout_point_Ex(f_point_Ex);
	ofstream fout_point_Ey(f_point_Ey);
	ofstream fout_point_Ez(f_point_Ez);

	ofstream fout_line_E(f_line_E);

	ofstream fout_Yface_Ex(f_Yface_Ex);
	ofstream fout_Yface_Ey(f_Yface_Ey);
	ofstream fout_Yface_Ez(f_Yface_Ez);
	//ofstream fout_Yface_Et(f_Yface_Et);

	ofstream fout_Zface_Ex(f_Zface_Ex);
	ofstream fout_Zface_Ey(f_Zface_Ey);
	ofstream fout_Zface_Ez(f_Zface_Ez);
	//ofstream fout_Zface_Et(f_Zface_Et);
	
	fout_point_Ex << scientific;
	fout_point_Ex.precision(4);
	fout_point_Ey << scientific;
	fout_point_Ey.precision(4);
	fout_point_Ez << scientific;
	fout_point_Ez.precision(4);

	fout_line_E << scientific;
	fout_line_E.precision(4);

	fout_Yface_Ex << scientific;
	fout_Yface_Ex.precision(4);
	fout_Yface_Ey << scientific;
	fout_Yface_Ey.precision(4);
	fout_Yface_Ez << scientific;
	fout_Yface_Ez.precision(4);

	fout_Zface_Ex << scientific;
	fout_Zface_Ex.precision(4);
	fout_Zface_Ey << scientific;
	fout_Zface_Ey.precision(4);
	fout_Zface_Ez << scientific;
	fout_Zface_Ez.precision(4);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//BEGIN TIME STEP
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::cout << "210527begin time-stepping" << std::endl;
	//添加时间循环
	for (int n = 0; n <= tmax; n++) {
		std::cout << "time n=" << n << std::endl;
		//t=(n-1/2)dt
		Einc.update1D_Hinc();
		// t=n*dt
		Einc.update1D_Einc(n);

		//*********************************************************************
		//******************         E t=n*dt              ********************
		update3D_E();   //t=n*dt
		update3D_CPML_psiE();
		Einc.add_TFSF_Box_E_CPML(ce_x_1, ce_x_2, ce_y_1, ce_y_2, ce_z_1, ce_z_2,
			psi_Exy_1, psi_Exy_2, psi_Exz_1, psi_Exz_2, psi_Eyx_1, psi_Eyx_2,
			psi_Eyz_1, psi_Eyz_2, psi_Ezx_1, psi_Ezx_2, psi_Ezy_1, psi_Ezy_2);
		update3D_CPML_E();
		Einc.add_TFSF_Box_E(CB, ob, den_ex, den_ey, den_ez, Ex, Ey, Ez);
		//Einc.add_TSFS_Z2_E(dx, dy, dz, CB, ob, den_ez, Ex, Ey);
		
		//*********************************************************************
		//******************       H t=(n+0.5)*dt           ********************
		update3D_H();
		update3D_CPML_psiH();
		Einc.add_TFSF_Box_H_CPML(ch_x_1, ch_x_2, ch_y_1, ch_y_2, ch_z_1, ch_z_2,
			psi_Hxy_1, psi_Hxy_2, psi_Hxz_1, psi_Hxz_2, psi_Hyx_1, psi_Hyx_2,
			psi_Hyz_1, psi_Hyz_2, psi_Hzx_1, psi_Hzx_2, psi_Hzy_1, psi_Hzy_2);
		update3D_CPML_H();	
		Einc.add_TFSF_Box_H(den_hx, den_hy, den_hz, Hx, Hy, Hz);
		//Einc.add_TSFS_Z2_H(dx, dy, dz, den_hz, Hx, Hy);


		//******************************************************************************
		//  需要修改
		if (n >= (tmax - 40)) {
			for (int i = Imin; i <= Imax-1; i++) {
				for (int j = Jmin; j <= Jmax-1; j++) {
					for (int k = Kmin; k <= Kmax-1; k++) {
						//Ex幅值
						temp = abs(Ex(i, j, k));
						if (temp > Emx(i, j,k)) {
							Emx(i, j, k) = temp;
						}
						temp = abs(Ey(i, j, k));
						if (temp > Emy(i, j, k)) {
							Emy(i, j, k) = temp;
						}
						//Ez幅值
						temp = abs(Ez(i, j, k));
						if (temp > Emz(i, j, k)) {
							Emz(i, j, k) = temp;
						}
					}//k
				}//j
			}//i
		}//n

		fout_point_Ex << setw(11) << n*dt*1.0E6
			<< setw(14) << Ex(0, 0, 0) << setw(14) << Ex(5, 5, 5) << setw(14) << Ex(0, 0, 10) << endl;
		fout_point_Ey << setw(11) << n*dt*1.0E6
			<< setw(14) << Ey(0, 0, 0) << setw(14) << Ey(5, 5, 5) << setw(14) << Ey(0, 0, 10) << endl;
		fout_point_Ez << setw(11) << n*dt*1.0E6
			<< setw(14) << Ez(0, 0, 0) << setw(14) << Ez(5, 5, 5) << setw(14) << Ez(0, 0, 10) << endl;
	}  // loop time
	
	
	////y面  j=0  Ex
	for (int k = KsMin; k <= KsMax -1; k++) {
		for (int i = IsMin; i <= IsMax - 1; i++) {
			tempx = (Emx(i, 0, k) + Emx(i, 1, k) + Emx(i, 0, k + 1) + Emx(i, 1, k + 1)) / 4.0;
			tempy = (Emy(i, 0, k) + Emy(i + 1, 0, k) + Emy(i, 0, k + 1) + Emy(i + 1, 0, k + 1)) / 4.0;
			tempz = (Emz(i, 0, k) + Emz(i+1 , 0, k) + Emz(i, 1, k) + Emz(i + 1, 1, k)) / 4.0;
			//temp = sqrt(tempx*tempx + tempy*tempy + tempz*tempz);

			fout_Yface_Ex << setw(14) << tempx;
			fout_Yface_Ey << setw(14) << tempy;
			fout_Yface_Ez << setw(14) << tempz;
			//fout_Yface_Et << setw(14) << temp;
		}
		fout_Yface_Ex << endl;
		fout_Yface_Ey << endl;
		fout_Yface_Ez << endl;
		//fout_Yface_Et << endl;
	}

	////z面 Ex
	int knum = 10;   //z方向观测面
	for (int j = JsMin; j <= JsMax -1; j++) {
		for (int i = IsMin; i <= IsMax - 1; i++) {
			tempx = (Emx(i, j, knum) + Emx(i, j + 1, knum) + Emx(i, j, knum + 1) + Emx(i, j + 1, knum + 1)) / 4.0;
			tempy = (Emy(i, j, knum) + Emy(i + 1, j, knum) + Emy(i, j, knum + 1) + Emy(i + 1, j, knum + 1)) / 4.0;
			tempz = (Emz(i, j, knum) + Emz(i + 1, j, knum) + Emz(i, j + 1, knum) + Emz(i + 1, j + 1, knum)) / 4.0;
			//temp = sqrt(tempx*tempx + tempy*tempy + tempz*tempz);

			fout_Zface_Ex << setw(14) << tempx;
			fout_Zface_Ey << setw(14) << tempy;
			fout_Zface_Ez << setw(14) << tempz;
			//fout_Zface_Et << setw(14) << temp;
		}
		fout_Zface_Ex << endl;
		fout_Zface_Ey << endl;
		fout_Zface_Ez << endl;
		//fout_Zface_Et << endl;
	}

	////z方向i=0 j=0
	for (int k = KsMin; k <= KsMax; k++) {
		fout_line_E << setw(11) << k*dz << setw(14) << Emx(0, 0, k)
			<< setw(14) << Emy(0, 0, k) << setw(14) << Emz(0, 0, k) << endl;
	}

	fout_point_Ex.close();
	fout_point_Ey.close();
	fout_point_Ez.close();
	fout_line_E.close();
	fout_Yface_Ex.close();
	fout_Yface_Ey.close();
	fout_Yface_Ez.close();
	//fout_Yface_Et.close();
	fout_Zface_Ex.close();
	fout_Zface_Ey.close();
	fout_Zface_Ez.close();
	//fout_Zface_Et.close();
}


void FDTD::IsMedia()
{
	////////start 粗糙地面建模
	//std::string fname = "read_data/L=500_h=0.5m_l=15m_N=1000_Lm=200m-Hm=57.74.dat";
	//ifstream fin_sur(fname);

	//if (!fin_sur.is_open()) {
	//	cout << fname << ":  open fail" << endl;
	//}

	//Matrix<double>ff(IsMin, IsMax, JsMin, JsMax);
	//for (int j = JsMin + 1; j <= JsMax; j++) {
	//	for (int i = IsMin + 1; i <= IsMax; i++) {
	//		fin_sur >> ff(i, j);
	//	}
	//}

	//for (int i = IsMin + 1; i <= IsMax; i++) {
	//	ff(i, JsMin) = ff(i,JsMax);
	//}
	//for (int j = JsMin + 1; j <= JsMax; j++) {
	//	ff(IsMin, j) = ff(IsMax, j);
	//}
	//ff(IsMin, JsMin) = ff(IsMax,JsMax);


	//double zz, zave;
	//for (int i = Imin; i <= Imax - 1; i++) {
	//	for (int j = Jmin; j <= Jmax - 1; j++) {
	//		for (int k = Kmin; k <= Kmax - 1; k++) {
	//			zz = (k + 0.5)*dz;
	//			if (i >= IsMin && i < IsMax && j >= JsMin && j < JsMax) {
	//				zave = (ff(i, j) + ff(i + 1, j) + ff(i, j + 1) + ff(i + 1, j + 1)) / 4.0;
	//				if (zz < zave)
	//					ob(i, j, k) = 1;
	//				else
	//					ob(i, j, k) = 0;
	//			}
	//			else {
	//				if (zz < 0.0)
	//					ob(i, j, k) = 1;
	//				else
	//					ob(i, j, k) = 0;
	//			}
	//		} //loop k
	//	}//loop j
	//}//loop i

	//string str_sur = "zz_sur.dat";
	//ofstream fout(str_sur);
	//fout.setf(ios_base::fixed, ios_base::floatfield);
	//fout.precision(4);
	////*******************************************************
	//for (int j = JsMin; j <= JsMax; j++) {
	//	for (int i = IsMin; i <= IsMax; i++) {
	//		fout.width(7);
	//		fout << ff(i, j) << "   ";
	//	}
	//	fout << endl;
	//}
	////*******************************************************
	//fout.close();
	////////end 粗糙地面建模



	for (int i = Imin; i <= Imax - 1; i++) {
		for (int j = Jmin; j <= Jmax - 1; j++) {
			for (int k = Kmin; k <= Kmax - 1; k++) {
				if (k + 0.5 <= 0.0) {
					ob(i, j, k)= 1;
				}
				else {
					ob(i, j, k) = 0;
				}
			}
		}
	}
}