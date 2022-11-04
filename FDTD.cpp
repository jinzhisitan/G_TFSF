#include <iostream>
#include"physical_constants.h"
#include"FDTD.h"

using namespace std;

FDTD::FDTD()
{
	std::cout << "class FDTD: void constructor\n";
}

FDTD::FDTD(int IsMin_, int IsMax_, int JsMin_, int JsMax_, int KsMin_, int KsMax_,
	int xPML1, int xPML2, int yPML1, int yPML2, int zPML1, int zPML2, int tmax_,
	double dx_, double dy_, double dz_, double dt_, vector<double> epsr_, vector<double> sig_,
	double alpha, double thi, double phi, int ItMin, int ItMax, int JtMin, int JtMax, int KtMin, int KtMax, int IncStart, int IncEnd)
	:Einc(dx_, dy_, dz_, alpha, thi, phi, dt_, ItMin, ItMax, JtMin, JtMax, KtMin, KtMax, IsMin_, IsMax_, JsMin_, JsMax_, KsMin_, KsMax_, IncStart, IncEnd),
	pml(xPML1, xPML2, yPML1, yPML2, zPML1, zPML2, IsMin_, IsMax_, JsMin_, JsMax_, KsMin_, KsMax_)
{
	IsMin = IsMin_; IsMax = IsMax_;
	JsMin = JsMin_; JsMax = JsMax_;
	KsMin = KsMin_; KsMax = KsMax_;

	Imin = IsMin_ - xPML1;
	Imax = IsMax_ + xPML2;
	Jmin = JsMin_ - yPML1;
	Jmax = JsMax_ + yPML2;
	Kmin = KsMin_ - zPML1;
	Kmax = KsMax_ + zPML2;

	tmax = tmax_;

	dx = dx_;
	dy = dy_;
	dz = dz_;
	dt = dt_;

	epsr= epsr_;
	sig = sig_;
}

void FDTD::setUp()
{
	initialize();
	initCoefficients();

	pml.initializeCPML();
	pml.initCoefficientsCPML(epsr, dx, dy, dz, dt, den_ex, den_ey, den_ez, den_hx, den_hy, den_hz);
	
	Einc.initializeTSFS();
	Einc.initCoefficientsTSFS();
}

void FDTD::initialize()
{
	Ex.Allocate(Imin, Imax - 1, Jmin, Jmax, Kmin, Kmax);
	Ey.Allocate(Imin, Imax, Jmin, Jmax - 1, Kmin, Kmax);
	Ez.Allocate(Imin, Imax, Jmin, Jmax, Kmin, Kmax - 1);
	Hx.Allocate(Imin, Imax, Jmin, Jmax - 1, Kmin, Kmax - 1);
	Hy.Allocate(Imin, Imax - 1, Jmin, Jmax, Kmin, Kmax - 1);
	Hz.Allocate(Imin, Imax - 1, Jmin, Jmax - 1, Kmin, Kmax);

	den_ex.Allocate(Imin + 1, Imax - 1), den_ey.Allocate(Jmin + 1, Jmax - 1), den_ez.Allocate(Kmin + 1, Kmax - 1);
	den_hx.Allocate(Imin, Imax - 1), den_hy.Allocate(Jmin, Jmax - 1), den_hz.Allocate(Kmin, Kmax - 1);

	ob.Allocate(Imin, Imax - 1, Jmin, Jmax - 1, Kmin, Kmax - 1);

	if (epsr.size() == sig.size()) {
		int MediaNo = int(epsr.size());
		CA.Allocate(0, MediaNo - 1, 0, MediaNo - 1, 0, MediaNo - 1, 0, MediaNo - 1);
		CB.Allocate(0, MediaNo - 1, 0, MediaNo - 1, 0, MediaNo - 1, 0, MediaNo - 1);
	}
	else {
		cout << "Error: Number of dielectric parameters" << endl;
	}
}

void FDTD::initCoefficients()
{
	double temp;

	//FILL IN UPDATING COEFFICIENTS
	DA = 1.0;
	DB = dt / mu0;
	int MediaNo = int(epsr.size());

	double epsr_eff, sig_eff;
	for (int i = 0; i < MediaNo; i++) {
		for (int j = 0; j < MediaNo; j++) {
			for (int k = 0; k < MediaNo; k++) {
				for (int l = 0; l < MediaNo; l++) {
					epsr_eff = 0.25*(epsr[i] + epsr[j] + epsr[k] + epsr[l]);    //4个网格介质的平均
					sig_eff = 0.25*(sig[i] + sig[j] + sig[k] + sig[l]);
					temp = sig_eff*dt / (2.0*eps0*epsr_eff);
					CA(i, j, k, l) = (1 - temp) / (1 + temp);    //Eq(2-2-17)
					CB(i, j, k, l) = dt / eps0 / (epsr_eff + sig_eff*dt / (2.0*eps0)); //Eq(2-2-18)
				}
			}	
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//FILL IN DENOMINATORS FOR FIELD UPDATES
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (int i = Imin; i <= Imax-1; i++) {
		den_hx(i) = 1.0 / dx;
	}
	for (int j = Jmin; j <= Jmax-1; j++) {
		den_hy(j) = 1.0 / dy;
	}
	for (int k = Kmin; k <= Kmax-1; k++) {
		den_hz(k) = 1.0 / dz;
	}
	for (int i = Imin+1; i <= Imax-1; i++) {
		den_ex(i) = 1.0 / dx;
	}
	for (int j = Jmin+1; j <= Jmax-1; j++) {
		den_ey(j) = 1.0 / dy;
	}
	for (int k = Kmin+1; k <= Kmax-1; k++) {
		den_ez(k) = 1.0 / dz;
	}

	//计算cell介电常数
	IsMedia();

}

void FDTD::update3D_H()
{
	int i, j, k;
	//UPDATE Hx
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax; i++) {
			for (j = Jmin; j <= Jmax - 1; j++) {
				//mm = hxmed(i, j, k);  //计算
				Hx(i, j, k) = DA * Hx(i, j, k) + DB *
					((Ez(i, j, k) - Ez(i, j + 1, k))*den_hy(j) + (Ey(i, j, k + 1) - Ey(i, j, k))*den_hz(k));
			}
		}
	}   //k-loop

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//UPDATE Hy
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax - 1; i++) {
			for (j = Jmin; j <= Jmax; j++) {
				//mm = hymed(i, j, k);
				Hy(i, j, k) = DA * Hy(i, j, k) + DB * 
					((Ez(i + 1, j, k) - Ez(i, j, k))*den_hx(i) + (Ex(i, j, k) - Ex(i, j, k + 1))*den_hz(k));
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Hz
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax; k++) {
		for (i = Imin; i <= Imax - 1; i++) {
			for (j = Jmin; j <= Jmax - 1; j++) {
				//mm = hzmed(i, j, k);
				Hz(i, j, k) = DA * Hz(i, j, k) + DB * 
					((Ey(i, j, k) - Ey(i + 1, j, k))*den_hx(i) + (Ex(i, j + 1, k) - Ex(i, j, k))*den_hy(j));
			}
		}
	}

}

void FDTD::update3D_E()
{
	int i, j, k;
	double factor1, factor2;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ex
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin + 1; k <= Kmax - 1; k++) {
		for (i = Imin; i <= Imax - 1; i++) {
			for (j = Jmin + 1; j <= Jmax - 1; j++) {
				factor1 = CA(ob(i, j, k), ob(i, j, k - 1), ob(i, j - 1, k), ob(i, j - 1, k - 1));
				factor2 = CB(ob(i, j, k), ob(i, j, k - 1), ob(i, j - 1, k), ob(i, j - 1, k - 1));

				Ex(i, j, k) = factor1 * Ex(i, j, k) + factor2 * ((Hz(i, j, k) - Hz(i, j - 1, k))*den_ey(j) +
					(Hy(i, j, k - 1) - Hy(i, j, k))*den_ez(k));
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ey
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin + 1; k <= Kmax - 1; k++) {
		for (i = Imin + 1; i <= Imax - 1; i++) {
			for (j = Jmin; j <= Jmax - 1; j++) {
				factor1 = CA(ob(i, j, k),ob(i - 1, j, k),ob(i, j, k - 1),ob(i - 1, j, k - 1));
				factor2 = CB(ob(i, j, k),ob(i - 1, j, k),ob(i, j, k - 1),ob(i - 1, j, k - 1));

				Ey(i, j, k) = factor1 * Ey(i, j, k) + factor2 * ((Hz(i - 1, j, k) - Hz(i, j, k))*den_ex(i) +
					(Hx(i, j, k) - Hx(i, j, k - 1))*den_ez(k));
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	//UPDATE Ez
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for (k = Kmin; k <= Kmax - 1; k++) {
		for (i = Imin + 1; i <= Imax - 1; i++) {
			for (j = Jmin + 1; j <= Jmax - 1; j++) {
				factor1 = CA(ob(i, j, k),ob(i - 1, j, k),ob(i, j - 1, k),ob(i - 1, j - 1, k));
				factor2 = CB(ob(i, j, k),ob(i - 1, j, k),ob(i, j - 1, k),ob(i - 1, j - 1, k));

				Ez(i, j, k) = factor1 * Ez(i, j, k) + factor2 * ((Hy(i, j, k) - Hy(i - 1, j, k))*den_ex(i) +
					(Hx(i, j - 1, k) - Hx(i, j, k))*den_ey(j));
			}
		}
	}
}