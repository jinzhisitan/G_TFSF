#include "TFSF.h"
#include<iostream>
#include <fstream>
#include<string>

void TFSF::initCoefficientsTSFS()
{
	int temp;
	double attnEtemp, attnHtemp;


	std::string fname00 = "../read_data/f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18.dat";
	std::ifstream fin00(fname00);
	if (!fin00.is_open()) {
		std::cout << fname00 << ":  open fail" << std::endl;
	}

	std::string fname45 = "../read_data/thi=45_f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18OLD.dat";
	std::ifstream fin45(fname45);
	if (!fin45.is_open()) {
		std::cout << fname45 << ":  open fail" << std::endl;
	}

	AttnE00.clear();
	AttnH00.clear();
	while (fin00 >> temp && fin00 >> attnEtemp&& fin00 >> attnHtemp)
	{
		AttnE00.push_back(attnEtemp);
		AttnH00.push_back(attnHtemp);
	}

	AttnE45.clear();
	AttnH45.clear();
	while (fin45 >> temp && fin45 >> attnEtemp&& fin45 >> attnHtemp)
	{
		AttnE45.push_back(attnEtemp);
		AttnH45.push_back(attnHtemp);
	}


	fin00.close();
	fin45.close();
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
	//kx沿着x方向正方向。
	if (!bIsAddHalfX)
	{
		if (i < IsMin)
			AxPML = pow(AttnE45[IsMin - i], -1);
		else if (i > IsMax)
			AxPML = pow(AttnE45[i - IsMax], 1);
		else
			AxPML = 1.0;
	}
	else {
		if (i < IsMin)
			AxPML = pow(AttnH45[IsMin - i - 1], -1);
		else if (i >= IsMax)
			AxPML = pow(AttnH45[i - IsMax], 1);
		else
			AxPML = 1.0;
	}


	////////////////////////////////////////////
	//******** z 方向
	//kz入射方向总是向下。
	if (!bIsAddHalfZ)
	{
		if (k < KsMin)
			AzPML = pow(AttnE45[KsMin - k], 1);
		else if (k > KsMax)
			AzPML = pow(AttnE45[k - KsMax], -1);
		else
			AzPML = 1.0;
	}
	else
	{
		if (k < KsMin)
			AzPML = pow(AttnH45[KsMin - k - 1], 1);
		else if (k >= KsMax)
			AzPML = pow(AttnH45[k - KsMax], -1);
		else
			AzPML = 1.0;
	}



	////////////////////////////////////////////
	//******** y 方向
	//if (!bIsAddHalfY)
	//{
	//	if (j < JsMin)
	//		AyPML = pow(1.0, -k0y);
	//	else if (j > JsMax)
	//		AyPML = pow(1.0, k0y);
	//	else
	//		AyPML = 1.0;
	//}
	//else
	//{
	//	if (j < JsMin)
	//		AyPML = pow(1.0, -k0y);
	//	else if (j >= JsMax)
	//		AyPML = pow(1.0, k0y);
	//	else
	//		AyPML = 1.0;
	//}

	//入射面为XZ面。
	AyPML = 1.0;


	APML = AxPML*AyPML*AzPML;
	return APML;
}