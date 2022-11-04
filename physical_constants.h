#pragma once
#include<cmath>
//Input Fundamental Constants(MKS units)
const double pi = 3.14159265358979;
const double c0 = 2.99792458E8; // ����
const double mu0 = 4.0*pi*1.0E-7; //��մŵ���
const double eps0 = 1.0 / (c0*c0*mu0); //��ս�糣��
const double Z0 = sqrt(mu0 / eps0);  //��ղ��迹


const int MediaNo = 2;  //��ͬ���ʵ�������
const int NApml = 15;   //�ܳ��߽���������PML���������