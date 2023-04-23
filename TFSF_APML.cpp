#include "TFSF.h"
#include<iostream>
#include <fstream>
#include<string>
#include<cmath>
#include"physical_constants.h"

void TFSF::initCoefficientsTSFS()
{
	int temp;
	double attnEtemp, attnHtemp;
	double phaseEtemp, phaseHtemp;

	std::string fname00 = "../read_data/f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18.dat";
	std::ifstream fin00(fname00);
	if (!fin00.is_open()) {
		std::cout << fname00 << ":  open fail" << std::endl;
	}

	std::string fname30 = "../read_data/thi=30_f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18.dat";
	std::ifstream fin30(fname30);
	if (!fin30.is_open()) {
		std::cout << fname30 << ":  open fail" << std::endl;
	}

	std::string fname45 = "../read_data/thi=45_f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18.dat";
	std::ifstream fin45(fname45);
	if (!fin45.is_open()) {
		std::cout << fname45 << ":  open fail" << std::endl;
	}

	std::string fname60 = "../read_data/thi=60_f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18.dat";
	std::ifstream fin60(fname60);
	if (!fin60.is_open()) {
		std::cout << fname60 << ":  open fail" << std::endl;
	}

	std::string fnamePhase00 = "../read_data/f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18_phase.dat";
	std::ifstream finPhase00(fnamePhase00);
	if (!finPhase00.is_open()) {
		std::cout << fnamePhase00 << ":  open fail" << std::endl;
	}

	std::string fnamePhase30 = "../read_data/thi=30_f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18_phase.dat";
	std::ifstream finPhase30(fnamePhase30);
	if (!finPhase30.is_open()) {
		std::cout << fnamePhase30 << ":  open fail" << std::endl;
	}

	std::string fnamePhase45 = "../read_data/thi=45_f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18_phase.dat";
	std::ifstream finPhase45(fnamePhase45);
	if (!finPhase45.is_open()) {
		std::cout << fnamePhase45 << ":  open fail" << std::endl;
	}

	std::string fnamePhase60 = "../read_data/thi=60_f=300MHz_dx=2.5cm_m=3_a=0_k=1_N=20inv18_phase.dat";
	std::ifstream finPhase60(fnamePhase60);
	if (!finPhase45.is_open()) {
		std::cout << fnamePhase60 << ":  open fail" << std::endl;
	}

	AttnE00.clear();
	AttnH00.clear();
	AttnE90.clear();
	AttnH90.clear();
	while (fin00 >> temp && fin00 >> attnEtemp&& fin00 >> attnHtemp)
	{
		AttnE00.push_back(attnEtemp);
		AttnH00.push_back(attnHtemp);
		AttnE90.push_back(1.0);
		AttnH90.push_back(1.0);
	}

	AttnE30.clear();
	AttnH30.clear();
	while (fin30 >> temp && fin30 >> attnEtemp&& fin30 >> attnHtemp)
	{
		AttnE30.push_back(attnEtemp);
		AttnH30.push_back(attnHtemp);
	}

	AttnE45.clear();
	AttnH45.clear();
	while (fin45 >> temp && fin45 >> attnEtemp&& fin45 >> attnHtemp)
	{
		AttnE45.push_back(attnEtemp);
		AttnH45.push_back(attnHtemp);
	}

	AttnE60.clear();
	AttnH60.clear();
	while (fin60 >> temp && fin60 >> attnEtemp&& fin60 >> attnHtemp)
	{
		AttnE60.push_back(attnEtemp);
		AttnH60.push_back(attnHtemp);
	}

	phaseE00.clear();
	phaseH00.clear();
	phaseE90.clear();
	phaseH90.clear();
	while (finPhase00 >> temp && finPhase00 >> phaseEtemp&& finPhase00 >> phaseHtemp)
	{
		phaseE00.push_back(phaseEtemp);
		phaseH00.push_back(phaseHtemp);
		phaseE90.push_back(0.0);
		phaseH90.push_back(0.0);
	}

	phaseE30.clear();
	phaseH30.clear();
	while (finPhase30 >> temp && finPhase30 >> phaseEtemp&& finPhase30 >> phaseHtemp)
	{
		phaseE30.push_back(phaseEtemp);
		phaseH30.push_back(phaseHtemp);
	}


	phaseE45.clear();
	phaseH45.clear();
	while (finPhase45 >> temp && finPhase45 >> phaseEtemp&& finPhase45 >> phaseHtemp)
	{
		phaseE45.push_back(phaseEtemp);
		phaseH45.push_back(phaseHtemp);
	}

	phaseE60.clear();
	phaseH60.clear();
	while (finPhase60 >> temp && finPhase60 >> phaseEtemp&& finPhase60 >> phaseHtemp)
	{
		phaseE60.push_back(phaseEtemp);
		phaseH60.push_back(phaseHtemp);
	}


	fin00.close();
	fin30.close();
	fin45.close();
	fin60.close();
	finPhase00.close();
	finPhase30.close();
	finPhase45.close();
	finPhase60.close();
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
	
	//////////////////////////////////////////////
	//////   解析计算
	//////////////////////////////////////////////
	//double k0x = sin(thi)*cos(phi);
	//double k0y = sin(thi)*sin(phi);
	//double k0z = cos(thi);
	////******** x 方向
	//if (!bIsAddHalfX)
	//{
	//	if (i < IsMin)
	//		AxPML = pow(AttnE00[IsMin - i], -k0x);
	//	else if (i > IsMax)
	//		AxPML = pow(AttnE00[i - IsMax], k0x);
	//	else
	//		AxPML = 1.0;
	//}
	//else {
	//	if (i < IsMin)
	//		AxPML = pow(AttnH00[IsMin - i - 1], -k0x);
	//	else if (i >= IsMax)
	//		AxPML = pow(AttnH00[i - IsMax], k0x);
	//	else
	//		AxPML = 1.0;
	//}
	//////////////////////////////////////////////
	////******** y 方向
	//if (!bIsAddHalfY)
	//{
	//	if (j < JsMin)
	//		AyPML = pow(AttnE00[JsMin - j], -k0y);
	//	else if (j > JsMax)
	//		AyPML = pow(AttnE00[j - JsMax], k0y);
	//	else
	//		AyPML = 1.0;
	//}
	//else
	//{
	//	if (j < JsMin)
	//		AyPML = pow(AttnH00[JsMin - j - 1], -k0y);
	//	else if (j >= JsMax)
	//		AyPML = pow(AttnH00[j - JsMax], k0y);
	//	else
	//		AyPML = 1.0;
	//}
	//////////////////////////////////////////////
	////******** z 方向
	//if (!bIsAddHalfZ)
	//{
	//	if (k < KsMin)
	//		AzPML = pow(AttnE00[KsMin - k], -k0z);
	//	else if (k > KsMax)
	//		AzPML = pow(AttnE00[k - KsMax], k0z);
	//	else
	//		AzPML = 1.0;
	//}
	//else
	//{
	//	if (k < KsMin)
	//		AzPML = pow(AttnH00[KsMin - k - 1], -k0z);
	//	else if (k >= KsMax)
	//		AzPML = pow(AttnH00[k - KsMax], k0z);
	//	else
	//		AzPML = 1.0;
	//}
	//APML = AxPML*AyPML*AzPML;
	//return APML;


	//////////////////////////////////////////////
	//   数值计算
	//////////////////////////////////////////////
	std::vector<double> AttnEx = AttnE90;
	std::vector<double> AttnHx = AttnH90;
	std::vector<double> AttnEz = AttnE00;
	std::vector<double> AttnHz = AttnH00;

	//******** x 方向
	//kx沿着x方向正方向。
	if (!bIsAddHalfX)
	{
		if (i < IsMin)
			AxPML = pow(AttnEx[IsMin - i], -1);
		else if (i > IsMax)
			AxPML = pow(AttnEx[i - IsMax], 1);
		else
			AxPML = 1.0;
	}
	else {
		if (i < IsMin)
			AxPML = pow(AttnHx[IsMin - i - 1], -1);
		else if (i >= IsMax)
			AxPML = pow(AttnHx[i - IsMax], 1);
		else
			AxPML = 1.0;
	}
	////////////////////////////////////////////
	//******** z 方向
	//kz入射方向总是向下。
	if (!bIsAddHalfZ)
	{
		if (k < KsMin)
			AzPML = pow(AttnEz[KsMin - k], 1);
		else if (k > KsMax)
			AzPML = pow(AttnEz[k - KsMax], -1);
		else
			AzPML = 1.0;
	}
	else
	{
		if (k < KsMin)
			AzPML = pow(AttnHz[KsMin - k - 1], 1);
		else if (k >= KsMax)
			AzPML = pow(AttnHz[k - KsMax], -1);
		else
			AzPML = 1.0;
	}
	////////////////////////////////////////////
	//******** y 方向
	////入射面为XZ面。
	AyPML = 1.0;

	APML = AxPML*AyPML*AzPML;
	return APML;
}


double TFSF::phaseCorrection(double x, double y, double z)
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


	double phaseX, phaseY, phaseZ;

	//////////////////////////////////////////////
	//   解析计算
	//////////////////////////////////////////////
	//double k0x = sin(thi)*cos(phi);
	//double k0y = sin(thi)*sin(phi);
	//double k0z = cos(thi);
	////******** x 方向
	//if (!bIsAddHalfX)
	//{
	//	if (i < IsMin)
	//		phaseX = phaseE00[IsMin - i] * (-k0x);  
	//	else if (i > IsMax)
	//		phaseX = phaseE00[i - IsMax] * k0x;
	//	else
	//		phaseX = 0.0;
	//}
	//else {
	//	if (i < IsMin)
	//		phaseX = phaseH00[IsMin - i - 1] *(-k0x);
	//	else if (i >= IsMax)
	//		phaseX = phaseH00[i - IsMax] * k0x;
	//	else
	//		phaseX = 0.0;
	//}
	//////////////////////////////////////////////
	////******** y 方向
	//if (!bIsAddHalfY)
	//{
	//	if (j < JsMin)
	//		phaseY = phaseE00[JsMin - j] *(-k0y);
	//	else if (j > JsMax)
	//		phaseY = phaseE00[j - JsMax] * k0y;
	//	else
	//		phaseY = 0.0;
	//}
	//else
	//{
	//	if (j < JsMin)
	//		phaseY = phaseH00[JsMin - j - 1] * (-k0y);
	//	else if (j >= JsMax)
	//		phaseY = phaseH00[j - JsMax] * k0y;
	//	else
	//		phaseY = 0.0;
	//}
	//////////////////////////////////////////////
	////******** z 方向
	//if (!bIsAddHalfZ)
	//{
	//	if (k < KsMin)
	//		phaseZ = phaseE00[KsMin - k]* (-k0z);
	//	else if (k > KsMax)
	//		phaseZ = phaseE00[k - KsMax]* k0z;
	//	else
	//		phaseZ = 0.0;
	//}
	//else
	//{
	//	if (k < KsMin)
	//		phaseZ = phaseH00[KsMin - k - 1] * (-k0z);
	//	else if (k >= KsMax)
	//		phaseZ = phaseH00[k - KsMax] * k0z;
	//	else
	//		phaseZ = 0.0;
	//}
	//double phase = phaseX+ phaseY+ phaseZ;
	////phase = 0.0;
	//return phase;


	////////////////////////////////////////////////
	////   数值计算
	////////////////////////////////////////////////
	std::vector<double> phaseEx = phaseE90;
	std::vector<double> phaseHx = phaseH90;
	std::vector<double> phaseEz = phaseE00;
	std::vector<double> phaseHz = phaseH00;

	//******** x 方向
	if (!bIsAddHalfX)
	{
		if (i < IsMin)
			phaseX = phaseEx[IsMin - i] * (-1);
		else if (i > IsMax)
			phaseX = phaseEx[i - IsMax] * 1;
		else
			phaseX = 0.0;
	}
	else {
		if (i < IsMin)
			phaseX = phaseHx[IsMin - i - 1] * (-1);
		else if (i >= IsMax)
			phaseX = phaseHx[i - IsMax] * 1;
		else
			phaseX = 0.0;
	}
	////////////////////////////////////////////
	//******** z 方向
	if (!bIsAddHalfZ)
	{
		if (k < KsMin)
			phaseZ = phaseEz[KsMin - k] * 1;
		else if (k > KsMax)
			phaseZ = phaseEz[k - KsMax] * (-1);
		else
			phaseZ = 0.0;
	}
	else
	{
		if (k < KsMin)
			phaseZ = phaseHz[KsMin - k - 1] * 1;
		else if (k >= KsMax)
			phaseZ = phaseHz[k - KsMax] * (-1);
		else
			phaseZ = 0.0;
	}
	////////////////////////////////////////////
	//******** y 方向
	////入射面为XZ面。
	phaseY = 0.0;
	double phase = phaseX + phaseY + phaseZ;
	return phase;
}