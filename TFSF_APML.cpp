#include<string>
#include"physical_constants.h"
#include"TFSF.h"

//程序运行前必须提供吸收层衰减系数文件。
//fname需要更改为所需要的文件。
void TFSF::APML()
{
	std::string fname = "read_data/f=300MHz_dx=05cm_m=4_a=0_k=0_N=16.dat";
	ifstream fin(fname);
	
	if (!fin.is_open()) {
		cout << fname << ":  open fail" << endl;
	}


	int temp;
	for (int i = 0; i < NApml; i++) {
		fin >> temp >> AttnE[i] >> AttnH[i];
	}
	fin.close();
}