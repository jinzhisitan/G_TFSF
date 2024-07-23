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

	bAnalysis = false;
	CalAttnEH();
	CalPhaseEH();
}


void TFSF::CalAttnEH()
{
	AttnEx.resize(AttnE00.size());
	AttnHx.resize(AttnH00.size());
	AttnEy.resize(AttnE00.size());
	AttnHy.resize(AttnH00.size());
	AttnEz.resize(AttnE00.size());
	AttnHz.resize(AttnH00.size());
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double k0xAbs = fabs(k0x);
	double k0yAbs = fabs(k0y);
	double k0zAbs = fabs(k0z);

	if (bAnalysis)
	{//解析计算
		for (int i = 0; i < AttnE00.size(); i++) {
			AttnEx[i] = pow(AttnE00[i], k0xAbs);
			AttnHx[i] = pow(AttnH00[i], k0xAbs);
			AttnEy[i] = pow(AttnE00[i], k0yAbs);
			AttnHy[i] = pow(AttnH00[i], k0yAbs);
			AttnEz[i] = pow(AttnE00[i], k0zAbs);
			AttnHz[i] = pow(AttnH00[i], k0zAbs);
		}
	}
	else  //数值计算
	{
		double Cos30 = sqrt(3.0) / 2.0;
		double Cos45 = sqrt(2.0) / 2.0;
		double Cos60 = 0.5;
		double w1, w2;
		//////////////////////////////////////////////////////////////////////
		////   x方向
		//////////////////////////////////////////////////////////////////////
		if (k0xAbs >= Cos30) {
			w1 = (1.0 - k0xAbs) / (1.0 - Cos30);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEx[i] = pow(AttnE00[i], w2)*pow(AttnE30[i], w1);
				AttnHx[i] = pow(AttnH00[i], w2)*pow(AttnH30[i], w1);
			}
		}
		else if (k0xAbs >= Cos45) {
			w1 = (Cos30 - k0xAbs) / (Cos30 - Cos45);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEx[i] = pow(AttnE30[i], w2)*pow(AttnE45[i], w1);
				AttnHx[i] = pow(AttnH30[i], w2)*pow(AttnH45[i], w1);
			}
		}
		else if (k0xAbs >= Cos60) {
			w1 = (Cos45 - k0xAbs) / (Cos45 - Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEx[i] = pow(AttnE45[i], w2)*pow(AttnE60[i], w1);
				AttnHx[i] = pow(AttnH45[i], w2)*pow(AttnH60[i], w1);
			}
		}
		else {
			w1 = (Cos60 - k0xAbs) / (Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEx[i] = pow(AttnE60[i], w2)*pow(AttnE90[i], w1);
				AttnHx[i] = pow(AttnH60[i], w2)*pow(AttnH90[i], w1);
			}
		}
		//////////////////////////////////////////////////////////////////////
		////   y方向
		//////////////////////////////////////////////////////////////////////
		if (k0yAbs >= Cos30) {
			w1 = (1.0 - k0yAbs) / (1.0 - Cos30);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEy[i] = pow(AttnE00[i], w2)*pow(AttnE30[i], w1);
				AttnHy[i] = pow(AttnH00[i], w2)*pow(AttnH30[i], w1);
			}
		}
		else if (k0yAbs >= Cos45) {
			w1 = (Cos30 - k0yAbs) / (Cos30 - Cos45);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEy[i] = pow(AttnE30[i], w2)*pow(AttnE45[i], w1);
				AttnHy[i] = pow(AttnH30[i], w2)*pow(AttnH45[i], w1);
			}
		}
		else if (k0yAbs >= Cos60) {
			w1 = (Cos45 - k0yAbs) / (Cos45 - Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEy[i] = pow(AttnE45[i], w2)*pow(AttnE60[i], w1);
				AttnHy[i] = pow(AttnH45[i], w2)*pow(AttnH60[i], w1);
			}
		}
		else {
			w1 = (Cos60 - k0yAbs) / (Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEy[i] = pow(AttnE60[i], w2)*pow(AttnE90[i], w1);
				AttnHy[i] = pow(AttnH60[i], w2)*pow(AttnH90[i], w1);
			}
		}
		//////////////////////////////////////////////////////////////////////
		////   z方向
		//////////////////////////////////////////////////////////////////////
		if (k0zAbs >= Cos30) {
			w1 = (1.0 - k0zAbs) / (1.0 - Cos30);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEz[i] = pow(AttnE00[i], w2)*pow(AttnE30[i], w1);
				AttnHz[i] = pow(AttnH00[i], w2)*pow(AttnH30[i], w1);
			}
		}
		else if (k0zAbs >= Cos45) {
			w1 = (Cos30 - k0zAbs) / (Cos30 - Cos45);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEz[i] = pow(AttnE30[i], w2)*pow(AttnE45[i], w1);
				AttnHz[i] = pow(AttnH30[i], w2)*pow(AttnH45[i], w1);
			}
		}
		else if (k0zAbs >= Cos60) {
			w1 = (Cos45 - k0zAbs) / (Cos45 - Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEz[i] = pow(AttnE45[i], w2)*pow(AttnE60[i], w1);
				AttnHz[i] = pow(AttnH45[i], w2)*pow(AttnH60[i], w1);
			}
		}
		else {
			w1 = (Cos60 - k0zAbs) / (Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < AttnE00.size(); i++) {
				AttnEz[i] = pow(AttnE60[i], w2)*pow(AttnE90[i], w1);
				AttnHz[i] = pow(AttnH60[i], w2)*pow(AttnH90[i], w1);
			}
		}
	}
}

void TFSF::CalPhaseEH()
{
	phaseEx.resize(phaseE00.size());
	phaseHx.resize(phaseH00.size());
	phaseEy.resize(phaseE00.size());
	phaseHy.resize(phaseH00.size());
	phaseEz.resize(phaseE00.size());
	phaseHz.resize(phaseH00.size());
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double k0xAbs = fabs(k0x);
	double k0yAbs = fabs(k0y);
	double k0zAbs = fabs(k0z);

	if (bAnalysis)
	{//解析计算
		for (int i = 0; i < phaseE00.size(); i++) {
			phaseEx[i] = phaseE00[i] * k0xAbs;
			phaseHx[i] = phaseH00[i] * k0xAbs;
			phaseEy[i] = phaseE00[i] * k0yAbs;
			phaseHy[i] = phaseH00[i] * k0yAbs;
			phaseEz[i] = phaseE00[i] * k0zAbs;
			phaseHz[i] = phaseH00[i] * k0zAbs;
		}
	}
	else  //数值计算
	{
		double Cos30 = sqrt(3.0) / 2.0;
		double Cos45 = sqrt(2.0) / 2.0;
		double Cos60 = 0.5;
		double w1, w2;
		//////////////////////////////////////////////////////////////////////
		////   x方向
		//////////////////////////////////////////////////////////////////////
		if (k0xAbs >= Cos30) {
			w1 = (1.0 - k0xAbs) / (1.0 - Cos30);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEx[i] = (phaseE00[i] * w2) + (phaseE30[i] * w1);
				phaseHx[i] = (phaseH00[i] * w2) + (phaseH30[i] * w1);
			}
		}
		else if (k0xAbs >= Cos45) {
			w1 = (Cos30 - k0xAbs) / (Cos30 - Cos45);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEx[i] = (phaseE30[i], w2) + (phaseE45[i], w1);
				phaseHx[i] = (phaseH30[i], w2) + (phaseH45[i], w1);
			}
		}
		else if (k0xAbs >= Cos60) {
			w1 = (Cos45 - k0xAbs) / (Cos45 - Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEx[i] = (phaseE45[i], w2) + (phaseE60[i], w1);
				phaseHx[i] = (phaseH45[i], w2) + (phaseH60[i], w1);
			}
		}
		else {
			w1 = (Cos60 - k0xAbs) / (Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEx[i] = (phaseE60[i], w2) + (phaseE90[i], w1);
				phaseHx[i] = (phaseH60[i], w2) + (phaseH90[i], w1);
			}
		}
		//////////////////////////////////////////////////////////////////////
		////   y方向
		//////////////////////////////////////////////////////////////////////
		if (k0yAbs >= Cos30) {
			w1 = (1.0 - k0yAbs) / (1.0 - Cos30);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEy[i] = (phaseE00[i], w2) + (phaseE30[i], w1);
				phaseHy[i] = (phaseH00[i], w2) + (phaseH30[i], w1);
			}
		}
		else if (k0yAbs >= Cos45) {
			w1 = (Cos30 - k0yAbs) / (Cos30 - Cos45);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEy[i] = (phaseE30[i], w2) + (phaseE45[i], w1);
				phaseHy[i] = (phaseH30[i], w2) + (phaseH45[i], w1);
			}
		}
		else if (k0yAbs >= Cos60) {
			w1 = (Cos45 - k0yAbs) / (Cos45 - Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEy[i] = (phaseE45[i], w2) + (phaseE60[i], w1);
				phaseHy[i] = (phaseH45[i], w2) + (phaseH60[i], w1);
			}
		}
		else {
			w1 = (Cos60 - k0yAbs) / (Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEy[i] = (phaseE60[i], w2) + (phaseE90[i], w1);
				phaseHy[i] = (phaseH60[i], w2) + (phaseH90[i], w1);
			}
		}
		//////////////////////////////////////////////////////////////////////
		////   z方向
		//////////////////////////////////////////////////////////////////////
		if (k0zAbs >= Cos30) {
			w1 = (1.0 - k0zAbs) / (1.0 - Cos30);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEz[i] = (phaseE00[i], w2) + (phaseE30[i], w1);
				phaseHz[i] = (phaseH00[i], w2) + (phaseH30[i], w1);
			}
		}
		else if (k0zAbs >= Cos45) {
			w1 = (Cos30 - k0zAbs) / (Cos30 - Cos45);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEz[i] = (phaseE30[i], w2) + (phaseE45[i], w1);
				phaseHz[i] = (phaseH30[i], w2) + (phaseH45[i], w1);
			}
		}
		else if (k0zAbs >= Cos60) {
			w1 = (Cos45 - k0zAbs) / (Cos45 - Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEz[i] = (phaseE45[i], w2) + (phaseE60[i], w1);
				phaseHz[i] = (phaseH45[i], w2) + (phaseH60[i], w1);
			}
		}
		else {
			w1 = (Cos60 - k0zAbs) / (Cos60);
			w2 = 1.0 - w1;
			for (int i = 0; i < phaseE00.size(); i++) {
				phaseEz[i] = (phaseE60[i], w2) + (phaseE90[i], w1);
				phaseHz[i] = (phaseH60[i], w2) + (phaseH90[i], w1);
			}
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
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double k0xSign = k0x > 0.0 ? 1.0 : -1.0;
	double k0ySign = k0y > 0.0 ? 1.0 : -1.0;
	double k0zSign = k0z > 0.0 ? 1.0 : -1.0;


	//******** x 方向
	if (!bIsAddHalfX)
	{
		if (i < IsMin)
			AxPML = pow(AttnEx[IsMin - i], -k0xSign);
		else if (i > IsMax)
			AxPML = pow(AttnEx[i - IsMax], k0xSign);
		else
			AxPML = 1.0;
	}
	else {
		if (i < IsMin)
			AxPML = pow(AttnHx[IsMin - i - 1], -k0xSign);
		else if (i >= IsMax)
			AxPML = pow(AttnHx[i - IsMax], k0xSign);
		else
			AxPML = 1.0;
	}
	////////////////////////////////////////////
	//******** y 方向
	if (!bIsAddHalfY)
	{
		if (j < JsMin)
			AyPML = pow(AttnEy[JsMin - j], -k0ySign);
		else if (j > JsMax)
			AyPML = pow(AttnEy[j - JsMax], k0ySign);
		else
			AyPML = 1.0;
	}
	else
	{
		if (j < JsMin)
			AyPML = pow(AttnHy[JsMin - j - 1], -k0ySign);
		else if (j >= JsMax)
			AyPML = pow(AttnHy[j - JsMax], k0ySign);
		else
			AyPML = 1.0;
	}
	////////////////////////////////////////////
	//******** z 方向
	if (!bIsAddHalfZ)
	{
		if (k < KsMin)
			AzPML = pow(AttnEz[KsMin - k], -k0zSign);
		else if (k > KsMax)
			AzPML = pow(AttnEz[k - KsMax], k0zSign);
		else
			AzPML = 1.0;
	}
	else
	{
		if (k < KsMin)
			AzPML = pow(AttnHz[KsMin - k - 1], -k0zSign);
		else if (k >= KsMax)
			AzPML = pow(AttnHz[k - KsMax], k0zSign);
		else
			AzPML = 1.0;
	}
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
	double k0x = sin(thi)*cos(phi);
	double k0y = sin(thi)*sin(phi);
	double k0z = cos(thi);
	double k0xSign = k0x > 0.0 ? 1.0 : -1.0;
	double k0ySign = k0y > 0.0 ? 1.0 : -1.0;
	double k0zSign = k0z > 0.0 ? 1.0 : -1.0;

	//******** x 方向
	if (!bIsAddHalfX)
	{
		if (i < IsMin)
			phaseX = phaseE00[IsMin - i] * (-k0xSign);
		else if (i > IsMax)
			phaseX = phaseE00[i - IsMax] * k0xSign;
		else
			phaseX = 0.0;
	}
	else {
		if (i < IsMin)
			phaseX = phaseH00[IsMin - i - 1] *(-k0xSign);
		else if (i >= IsMax)
			phaseX = phaseH00[i - IsMax] * k0xSign;
		else
			phaseX = 0.0;
	}
	////////////////////////////////////////////
	//******** y 方向
	if (!bIsAddHalfY)
	{
		if (j < JsMin)
			phaseY = phaseE00[JsMin - j] *(-k0ySign);
		else if (j > JsMax)
			phaseY = phaseE00[j - JsMax] * k0ySign;
		else
			phaseY = 0.0;
	}
	else
	{
		if (j < JsMin)
			phaseY = phaseH00[JsMin - j - 1] * (-k0ySign);
		else if (j >= JsMax)
			phaseY = phaseH00[j - JsMax] * k0ySign;
		else
			phaseY = 0.0;
	}
	////////////////////////////////////////////
	//******** z 方向
	if (!bIsAddHalfZ)
	{
		if (k < KsMin)
			phaseZ = phaseE00[KsMin - k]* (-k0zSign);
		else if (k > KsMax)
			phaseZ = phaseE00[k - KsMax]* k0zSign;
		else
			phaseZ = 0.0;
	}
	else
	{
		if (k < KsMin)
			phaseZ = phaseH00[KsMin - k - 1] * (-k0zSign);
		else if (k >= KsMax)
			phaseZ = phaseH00[k - KsMax] * k0zSign;
		else
			phaseZ = 0.0;
	}
	double phase = phaseX+ phaseY+ phaseZ;
	//phase = 0.0;
	return phase;
}