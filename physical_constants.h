#pragma once
#include<cmath>
//Input Fundamental Constants(MKS units)
const double pi = 3.14159265358979;
const double c0 = 2.99792458E8; // 光速
const double mu0 = 4.0*pi*1.0E-7; //真空磁导率
const double eps0 = 1.0 / (c0*c0*mu0); //真空介电常数
const double Z0 = sqrt(mu0 / eps0);  //真空波阻抗


const int MediaNo = 2;  //不同介质的数量。
const int NApml = 15;   //总场边界条件深入PML层的最大层数