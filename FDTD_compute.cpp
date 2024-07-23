#include <iostream>
#include <fstream>
#include<iomanip>
#include<cmath>
#include<string>
#include"physical_constants.h"
#include"FDTD.h"

using namespace std;

void FDTD::compute()
{
	double tempx, tempy, tempz, temp;
	Matrix<double> Emx(Imin, Imax - 1, Jmin, Jmax - 1, Kmin, Kmax - 1);
	Matrix<double> Emy(Imin, Imax - 1, Jmin, Jmax - 1, Kmin, Kmax - 1);
	Matrix<double> Emz(Imin, Imax - 1, Jmin, Jmax - 1, Kmin, Kmax - 1);

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
	std::cout << "210527FDTDbegin time-stepping" << std::endl;
	for (int n = 0; n <= tmax; n++) {
		std::cout << "time n=" << n << std::endl;
		//*********************************************************************
		//******************         Begin Main            ********************
		//*********************************************************************
		//***t=(n-1/2)dt
		Einc.update1D_Hinc();
		//***t=n*dt
		Einc.update1D_Einc();
		Einc.addSource(n);

		//******************         E t=n*dt              ********************
		update3D_E();   //t=n*dt
		pml.update3D_CPML_psiE(Hx, Hy, Hz);
		Einc.add_TFSF_Box_E_CPML(pml.ce_x_1, pml.ce_x_2, pml.ce_y_1, pml.ce_y_2, pml.ce_z_1, pml.ce_z_2,
			pml.psi_Exy_1, pml.psi_Exy_2, pml.psi_Exz_1, pml.psi_Exz_2, pml.psi_Eyx_1, pml.psi_Eyx_2,
			pml.psi_Eyz_1, pml.psi_Eyz_2, pml.psi_Ezx_1, pml.psi_Ezx_2, pml.psi_Ezy_1, pml.psi_Ezy_2);
		//Einc.add_TFSF_Box_E_CPML_analysis(n-0.5, pml.ce_x_1, pml.ce_x_2, pml.ce_y_1, pml.ce_y_2, pml.ce_z_1, pml.ce_z_2,
		//	pml.psi_Exy_1, pml.psi_Exy_2, pml.psi_Exz_1, pml.psi_Exz_2, pml.psi_Eyx_1, pml.psi_Eyx_2,
		//	pml.psi_Eyz_1, pml.psi_Eyz_2, pml.psi_Ezx_1, pml.psi_Ezx_2, pml.psi_Ezy_1, pml.psi_Ezy_2);
		pml.update3D_CPML_E(CB, ob, Ex, Ey, Ez);
		Einc.add_TFSF_Box_E(CB, ob, den_ex, den_ey, den_ez, Ex, Ey, Ez);
		//Einc.add_TSFS_Box_E_analysis(n-0.5, CB, ob, den_ex, den_ey, den_ez, Ex, Ey, Ez);
		//Einc.add_TFSF_Z2_E(Ex, Ey);

		//******************       H t=(n+0.5)*dt           ********************
		update3D_H();
		pml.update3D_CPML_psiH(Ex, Ey, Ez);
		Einc.add_TFSF_Box_H_CPML(pml.ch_x_1, pml.ch_x_2, pml.ch_y_1, pml.ch_y_2, pml.ch_z_1, pml.ch_z_2,
			pml.psi_Hxy_1, pml.psi_Hxy_2, pml.psi_Hxz_1, pml.psi_Hxz_2, pml.psi_Hyx_1, pml.psi_Hyx_2,
			pml.psi_Hyz_1, pml.psi_Hyz_2, pml.psi_Hzx_1, pml.psi_Hzx_2, pml.psi_Hzy_1, pml.psi_Hzy_2);
		//Einc.add_TFSF_Box_H_CPML_analysis(n, pml.ch_x_1, pml.ch_x_2, pml.ch_y_1, pml.ch_y_2, pml.ch_z_1, pml.ch_z_2,
		//	pml.psi_Hxy_1, pml.psi_Hxy_2, pml.psi_Hxz_1, pml.psi_Hxz_2, pml.psi_Hyx_1, pml.psi_Hyx_2,
		//	pml.psi_Hyz_1, pml.psi_Hyz_2, pml.psi_Hzx_1, pml.psi_Hzx_2, pml.psi_Hzy_1, pml.psi_Hzy_2);
		pml.update3D_CPML_H(DB, Hx, Hy, Hz);
		Einc.add_TFSF_Box_H(den_hx, den_hy, den_hz, Hx, Hy, Hz);
		//Einc.add_TSFS_Box_H_analysis(n, den_hx, den_hy, den_hz,Hx, Hy, Hz);
		//Einc.add_TFSF_Z2_H(Hx, Hy);

		//**********************************************************************
		//******************         End Main               ********************
		//*********************************************************************
		
		//double T = 2.0E-9;
		//Ez(0, 0, 0) = Ez(0, 0, 0) + dt / eps0 / dx / dx / dx*
		//	1.0E-10* (-2.0*  ((n*dt - 3.0*T) / T) / T)
		//	*exp(-pow((n*dt - 3.0*T) / T, 2)  );


		if (n >= (tmax - 50)) {
			for (int i = Imin; i <= Imax - 1; i++) {
				for (int j = Jmin; j <= Jmax - 1; j++) {
					for (int k = Kmin; k <= Kmax - 1; k++) {
						//Ex幅值
						temp = abs(Ex(i, j, k));
						if (temp > Emx(i, j, k)) {
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
			<< setw(14) << Ex(0, 0, 4)
			<< setw(14) << Ex(0, 0, 20)
			<< setw(14) << Ex(0, 0, 40)
			<< setw(14) << Ex(0, 0, 80)
			 << endl;
		fout_point_Ey << setw(11) << n*dt*1.0E6
			<< setw(14) << Ey(0, 0, 2)
			<< setw(14) << Ey(0, 0, 10)
			<< setw(14) << Ey(0, 0, 20)
			<< setw(14) << Ey(0, 0, 22)
			<< setw(14) << Ey(0, 0, 25)
			<< endl;
		fout_point_Ez << setw(11) << n*dt*1.0E6
			<< setw(14) << Ez(0, 0, 2)
			<< setw(14) << Ez(0, 0, 20)
			<< setw(14) << Ez(0, 0, 40)
			<< setw(14) << Ez(0, 0, 80)
			<< endl;


		////y面  j=0  Ex
		//if (n == 300) {
		//for (int k = KsMin; k <= KsMax - 1; k++) {
		//	for (int i = IsMin; i <= IsMax - 1; i++) {
		//		tempx = abs(Ex(i, 0, k));
		//		tempy = abs(Ey(i, 0, k));
		//		tempz = abs(Ez(i, 0, k));

		//		fout_Yface_Ex << setw(14) << tempx;
		//		fout_Yface_Ey << setw(14) << tempy;
		//		fout_Yface_Ez << setw(14) << tempz;
		//	}
		//	fout_Yface_Ex << endl;
		//	fout_Yface_Ey << endl;
		//	fout_Yface_Ez << endl;
		//}
		//}
	}  // loop time

	   ////y面  j=0  Ex
	for (int k = KsMin; k <= KsMax-5; k++) {
		for (int i = IsMin; i <= IsMax-1; i++) {
			tempx = Emx(i, 0, k);
			fout_Yface_Ex << setw(14) << tempx;
		}
		fout_Yface_Ex << endl;
	}


	for (int k = KsMin; k <= KsMax - 5; k++) {
		for (int i = IsMin; i <= IsMax; i++) {
			tempy = Emy(i, 0, k);
			fout_Yface_Ey << setw(14) << tempy;
		}
		fout_Yface_Ey << endl;
	}

	for (int k = KsMin; k <= KsMax - 5-1; k++) {
		for (int i = IsMin; i <= IsMax; i++) {
			tempz = Emz(i, 0, k);
			fout_Yface_Ez << setw(14) << tempz;
		}
		fout_Yface_Ez << endl;
	}



	////z面 Ex
	int knum = 10;   //z方向观测面
	for (int j = JsMin; j <= JsMax; j++) {
		for (int i = IsMin; i <= IsMax-1; i++) {
			tempx = Emx(i, j, knum);
			fout_Zface_Ex << setw(14) << tempx;
		}
		fout_Zface_Ex << endl;
	}

	for (int j = JsMin; j <= JsMax - 1; j++) {
		for (int i = IsMin; i <= IsMax; i++) {
			tempy = Emy(i, j, knum);
			fout_Zface_Ey << setw(14) << tempy;	
		}
		fout_Zface_Ey << endl;
	}


	for (int j = JsMin; j <= JsMax; j++) {
		for (int i = IsMin; i <= IsMax; i++) {
			tempz = Emz(i, j, knum);
			fout_Zface_Ez << setw(14) << tempz;
		}
		fout_Zface_Ez << endl;
	}


	//z方向i=0 j=0
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
	////  1  上层自由空间，下层地面。
	//for (int i = Imin; i <= Imax - 1; i++) {
	//	for (int j = Jmin; j <= Jmax - 1; j++) {
	//		for (int k = Kmin; k <= Kmax - 1; k++) {
	//			if (k + 0.5 <= 0.0) {
	//				ob(i, j, k) = 1;
	//			}
	//			else {
	//				ob(i, j, k) = 0;
	//			}
	//		}
	//	}
	//}


	////  2  上层自由空间，下层地面。
	for (int i = Imin; i <= Imax - 1; i++) {
		for (int j = Jmin; j <= Jmax - 1; j++) {
			for (int k = Kmin; k <= Kmax - 1; k++) {
				ob(i, j, k) = 0;
			}
		}
	}
	int ScaterGrid = 16;
	for (int i = IsMin - ScaterGrid; i <= IsMax + ScaterGrid - 1; i++) {
		for (int j = JsMin - ScaterGrid; j <= JsMax + ScaterGrid - 1; j++) {
			for (int k = KsMin - 17; k <= -1; k++) {
				//if ((k + 0.5) > -10.0) {    //第1层，厚度为
				//	ob(i, j, k) = 1;
				//}
				//else if ((k + 0.5) > -10.0){ //第2层，厚度为
				//	ob(i, j, k) = 2;
				//}
				//else if((k + 0.5) > -15.0){
				//	ob(i, j, k) = 3;
				//}
				//else if ((k + 0.5) > -20.0) {
				//	ob(i, j, k) = 4;
				//}
				//else {
				//	ob(i, j, k) = 2;
				//}	

	
				ob(i, j, k) = 1;
			}
		}
	}


	//int ScaterGrid = 16;
	//int IsMinTemp = IsMin - ScaterGrid, IsMaxTemp = IsMax + ScaterGrid,
	//	JsMinTemp = IsMinTemp, JsMaxTemp = IsMaxTemp;
	////首先全部设置为自由空间。
	//for (int i = Imin; i <= Imax - 1; i++) {
	//	for (int j = Jmin; j <= Jmax - 1; j++) {
	//		for (int k = Kmin; k <= Kmax - 1; k++) {
	//			ob(i, j, k) = 0;
	//		}
	//	}
	//}
	////////start 粗糙地面建模
	//std::string fname = "../read_data/L=10.8_h=0.1m_l=1m_N=432.dat";
	//ifstream fin_sur(fname);
	//if (!fin_sur.is_open()) {
	//	std::cout << fname << ":  open fail" << endl;
	//}
	//
	//Matrix<double>ff(IsMinTemp, IsMaxTemp, JsMinTemp, JsMaxTemp);
	//
	//for (int j = JsMinTemp + 1; j <= JsMaxTemp; j++) {
	//	for (int i = IsMinTemp + 1; i <= IsMaxTemp; i++) {
	//		fin_sur >> ff(i, j);
	//	}
	//}

	//for (int i = IsMinTemp + 1; i <= IsMaxTemp; i++) {
	//	ff(i, JsMinTemp) = ff(i, JsMaxTemp);
	//}
	//for (int j = JsMinTemp + 1; j <= JsMaxTemp; j++) {
	//	ff(IsMinTemp, j) = ff(IsMaxTemp, j);
	//}
	//ff(IsMinTemp, JsMinTemp) = ff(IsMaxTemp, JsMaxTemp);

	//double zz, zave;
	//for (int i = IsMinTemp; i <= IsMaxTemp - 1; i++) {
	//	for (int j = JsMinTemp; j <= JsMaxTemp - 1; j++) {
	//		for (int k = KsMin - 17; k <= KsMax-5; k++) {
	//			zz = (k + 0.5)*dz;
	//			zave = (ff(i, j) + ff(i + 1, j) + ff(i, j + 1) + ff(i + 1, j + 1)) / 4.0;
	//			if (zz < zave)
	//				ob(i, j, k) = 1;
	//		} //loop k
	//	}//loop j
	//}//loop i

	//string str_sur = "zz_sur.dat";
	//ofstream fout(str_sur);
	//fout.setf(ios_base::fixed, ios_base::floatfield);
	//fout.precision(4);
	////*******************************************************
	//for (int j = JsMinTemp; j <= JsMaxTemp; j++) {
	//	for (int i = IsMinTemp; i <= IsMaxTemp; i++) {
	//		fout.width(7);
	//		fout << ff(i, j) << "   ";
	//	}
	//	fout << endl;
	//}
	////*******************************************************
	//fout.close();
	////////end 粗糙地面建模
}