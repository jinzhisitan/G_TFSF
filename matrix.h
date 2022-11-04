#pragma once

//PURPOSE:
//	 动态分配1-4维中任意维度的数组，并且每一维度的上下限可以为负数。
//NOTES:
//	1.若对象调用默认构造函数Matrix()，则需要进一步调用相应维度的Allocate()函数为数组分配空间和赋值。
//  2.数组内存排列方式先排列最高维数据。
//	3.需要调用重构符号operator()访问数组中的某个值。
//Record of revisions:
//   1. 06 June 2020 by Yang chao, Email: yangchaophy@foxmail.com
//      Northwest Institute of Nuclear Technology, Xi'an; 710024, China
//	 2. 27 May 2021 by Yang chao, Email: yangchaophy@foxmail.com
//      Northwest Institute of Nuclear Technology, Xi'an; 710024, China

#include<iostream>
using namespace std;

template<class T>
class Matrix {
private:
	int amin, amax; // 第1维度的上限和下限,可为负数（amax>=amin）。
	int	bmin, bmax; // 第2维度的上限和下限,可为负数（bmax>=bmin）。
	int	cmin, cmax; // 第3维度的上限和下限,可为负数（cmax>=cmin）。
	int dmin, dmax; // 第4维度的上限和下限,可为负数（dmax>=dmin）。

	T* data;   //动态分配的数组
public:
	//默认构造函数。若调用默认构造函数，则需要进一步调用相应维度的Allocate()分配数组空间和赋值。
	Matrix();  
	//为1维度数组分配空间和赋值。
    //a1:下限，a2:上限; val: 数组的值。
	Matrix(int a1, int a2, T val=0);
	Matrix(int a1, int a2, int b1, int b2, T val = 0);
	Matrix(int a1, int a2, int b1, int b2, int c1, int c2, T val = 0);
	//为4维度数组分配空间和赋值。
    //a1:第1维下限，a2:第1维上限;  b1:第2维下限，b2:第2维上限;
	//c1:第3维下限，c2:第3维上限;  d1:第4维下限，d2:第4维上限;
	// val: 数组的值。
	Matrix(int a1, int a2, int b1, int b2, int c1, int c2, int d1, int d2, T val = 0);
	//copy constructor
	Matrix(const Matrix &orig);
	//重构对象赋值运算符。
	Matrix & operator=(const Matrix &orig);
	//deconstructor
	~Matrix();

	//1维数组分配空间。
	//a1:下限，a2:上限; val: 数组的值。
	void Allocate(int a1, int a2, T val=0);
	void Allocate(int a1, int a2, int b1, int b2, T val=0);
	void Allocate(int a1, int a2, int b1, int b2, int c1, int c2, T val=0);
	//为4维度数组分配空间和赋值。
	//a1:第1维下限，a2:第1维上限;  b1:第2维下限，b2:第2维上限;
	//c1:第3维下限，c2:第3维上限;  d1:第4维下限，d2:第4维上限;
	//val: 数组的值。
	void Allocate(int a1, int a2, int b1, int b2, int c1, int c2, int d1, int d2, T val=0);
	//initial array to @c initval
	//@param initval
	void initArray(T val = 0);

	//overload the parents。
	//返回1维数组中的data[a-amin]值
	T & operator()(int a);
	T & operator()(int a, int b);
	T & operator()(int a, int b, int c);
	//返回4维数组中的data[a-amin][b - amin][c - amin][d - amin]
	T & operator()(int a, int b, int c, int d);

	const T & operator()(int a)const;
	const T & operator()(int a, int b)const;
	const T & operator()(int a, int b, int c)const;
	const T & operator()(int a, int b, int c, int d)const;
};

//default donstructor
//set @c data  nullptr
//array length zero
template<class T>
Matrix<T>::Matrix()
{
	amin = 0; amax = 0;
	bmin = 0; bmax = 0;
	cmin = 0; cmax = 0;
	dmin = 0; dmax = 0;
	data = nullptr;
}

template<class T>
Matrix<T>::Matrix(int a1,int a2, T val)
{
	if (a2 >= a1) {
		amin = a1; amax = a2;
		bmin = 0; bmax = 0;
		cmin = 0; cmax = 0;
		dmin = 0; dmax = 0;

		int asize = a2 - a1 + 1;
		data = new T[asize];

		initArray(val);
	}
	else {
		std::cout << "Error in Matrix::Allocate(amin, amax) " << std::endl;
	}
}

template<class T>
Matrix<T>::Matrix(int a1, int a2, int b1, int b2, T val)
{
	if (a2 >= a1 && b2 >= b1) {
		amin = a1; amax = a2;
		bmin = b1; bmax = b2;
		cmin = 0; cmax = 0;
		dmin = 0; dmax = 0;

		int asize = a2 - a1 + 1;
		int bsize = b2 - b1 + 1;
		data = new T[asize * bsize];

		initArray(val);
	}
	else {
		std::cout << "Error in Matrix::Allocate(amin, amax, bmin, bmax) " << std::endl;
	}
}

template<class T>
Matrix<T>::Matrix(int a1, int a2, int b1, int b2, int c1, int c2, T val)
{
	if (a2 >= a1 && b2 >= b1 && c2 >= c1) {
		amin = a1; amax = a2;
		bmin = b1; bmax = b2;
		cmin = c1; cmax = c2;
		dmin = 0; dmax = 0;

		int asize = a2 - a1 + 1;
		int bsize = b2 - b1 + 1;
		int csize = c2 - c1 + 1;
		data = new T[asize * bsize * csize];

		initArray(val);
	}
	else {
		std::cout << "Error in Matrix::Allocate(amin, amax, bmin, bmax, cmin, cmax) " << std::endl;
	}
}

template<class T>
Matrix<T>::Matrix(int a1, int a2, int b1, int b2, int c1, int c2, int d1, int d2, T val)
{
	if (a2 >= a1 && b2 >= b1 && c2 >= c1 && d2 >= d1) {
		amin = a1; amax = a2;
		bmin = b1; bmax = b2;
		cmin = c1; cmax = c2;
		dmin = d1; dmax = d2;

		int asize = a2 - a1 + 1;
		int bsize = b2 - b1 + 1;
		int csize = c2 - c1 + 1;
		int dsize = d2 - d1 + 1;
		data = new T[asize * bsize * csize * dsize];

		initArray(val);
	}
	else {
		std::cout << "Error in Matrix::allocate(amin, amax, bmin, bmax, cmin, cmax, dmin, dmax,) " << std::endl;
	}
}

template<class T>
Matrix<T>::Matrix(const Matrix<T>&orig)
{
	amin = orig.amin; amax = orig.amax;
	bmin = orig.bmin; bmax = orig.bmax;
	cmin = orig.cmin; cmax = orig.cmax;
	dmin = orig.dmin; dmax = orig.dmax;

	int asize = orig.amax - orig.amin + 1;
	int bsize = orig.bmax - orig.bmin + 1;
	int csize = orig.cmax - orig.cmin + 1;
	int dsize = orig.dmax - orig.dmin + 1;
	data = new T[asize * bsize * csize * dsize];

	int n;
	for (int a = 0; a < asize; a++) {
		for (int b = 0; b < bsize; b++) {
			for (int c = 0; c < csize; c++) {
				for (int d = 0; d < dsize; d++) {
					n = a*bsize*csize*dsize + b*csize*dsize + c*dsize + d;
					data[n] = orig.data[n];
				}
			}
		}
	}
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &orig)
{
	if (this == &orig) {
		return *this;
	}
	delete[] data;

	amin = orig.amin; amax = orig.amax;
	bmin = orig.bmin; bmax = orig.bmax;
	cmin = orig.cmin; cmax = orig.cmax;
	dmin = orig.dmin; dmax = orig.dmax;

	int asize = orig.amax - orig.amin + 1;
	int bsize = orig.bmax - orig.bmin + 1;
	int csize = orig.cmax - orig.cmin + 1;
	int dsize = orig.dmax - orig.dmin + 1;
	data = new T[asize * bsize * csize * dsize];

	int n;
	for (int a = 0; a < asize; a++) {
		for (int b = 0; b < bsize; b++) {
			for (int c = 0; c < csize; c++) {
				for (int d = 0; d < dsize; d++) {
					n = a*bsize*csize*dsize + b*csize*dsize + c*dsize + d;
					data[n] = orig.data[n];
				}
			}
		}
	}
	return *this;
}

template<class T>
Matrix<T>::~Matrix()
{
	delete[]data;
}

template<class T>
void Matrix<T>::Allocate(int a1, int a2, T val)
{
	if (a2 >= a1) {
		if (data != nullptr) {
			delete[]data;
		}

		amin = a1; amax = a2;
		bmin = 0; bmax = 0;
		cmin = 0; cmax = 0;
		dmin = 0; dmax = 0;

		int asize = a2 - a1 + 1;
		data = new T[asize];
		
		initArray(val);
	}
	else {
		std::cout << "Error in Matrix::Allocate(amin, amax) " << std::endl;
	}
}

template<class T>
void Matrix<T>::Allocate(int a1, int a2, int b1, int b2, T val)
{
	if (a2 >= a1 && b2 >= b1) {
		if (data != nullptr) {
			delete[]data;
		}
		
		amin = a1; amax = a2;
		bmin = b1; bmax = b2;
		cmin = 0; cmax = 0;
		dmin = 0; dmax = 0;

		int asize = a2 - a1 + 1;
		int bsize = b2 - b1 + 1;
		data = new T[asize * bsize];
		
		initArray(val);
	}
	else {
		std::cout << "Error in Matrix::Allocate(amin, amax, bmin, bmax) " << std::endl;
	}
}

template<class T>
void Matrix<T>::Allocate(int a1, int a2, int b1, int b2, int c1, int c2, T val)
{
	if (a2 >= a1 && b2 >= b1 && c2 >= c1) {
		if (data != nullptr) {
			delete[]data;
		}

		amin = a1; amax = a2;
		bmin = b1; bmax = b2;
		cmin = c1; cmax = c2;
		dmin = 0; dmax = 0;

		int asize = a2 - a1 + 1;
		int bsize = b2 - b1 + 1;
		int csize = c2 - c1 + 1;
		data = new T[asize * bsize * csize];

		initArray(val);
	}
	else {
		std::cout << "Error in Matrix::Allocate(amin, amax, bmin, bmax, cmin, cmax) " << std::endl;
	}
}

template<class T>
void Matrix<T>::Allocate(int a1, int a2, int b1, int b2, int c1, int c2, int d1, int d2, T val)
{
	if (a2 >= a1 && b2 >= b1 && c2 >= c1 && d2>=d1) {
		if (data != nullptr) {
			delete[]data;
		}

		amin = a1; amax = a2;
		bmin = b1; bmax = b2;
		cmin = c1; cmax = c2;
		dmin = d1; dmax = d2;

		int asize = a2 - a1 + 1;
		int bsize = b2 - b1 + 1;
		int csize = c2 - c1 + 1;
		int dsize = d2 - d1 + 1;
		data = new T[asize * bsize * csize * dsize];
		
		initArray(val);
	}
	else {
		std::cout << "Error in Matrix::allocate(amin, amax, bmin, bmax, cmin, cmax, dmin, dmax,) " << std::endl;
	}
}

template<class T>
void Matrix<T>::initArray(T val)
{
	if (data == nullptr) {
		return;
	}
	
	int asize = amax - amin + 1;
	int bsize = bmax - bmin + 1;
	int csize = cmax - cmin + 1;
	int dsize = dmax - dmin + 1;

	for (int a = 0; a < asize; a++) {
		for (int b = 0; b < bsize; b++) {
			for (int c = 0; c < csize; c++) {
				for (int d = 0; d < dsize; d++) {
					data[a*bsize*csize*dsize + b*csize*dsize + c*dsize + d] = val;
				}
			}
		}
	}

}

template<class T>
T& Matrix<T>::operator()(int a)
{
	if (a >= amin&&a <= amax) {
		int n = a - amin;
		return data[n];
	}
	else {
		std::cout << "Error in index" << std::endl;
		return data[0];
	}
	
}

template<class T>
T& Matrix<T>::operator()(int a, int b)
{
	if (a >= amin&&a <= amax&&b >= bmin&&b <= bmax) {
		int bsize = bmax - bmin + 1;
		int n = (a - amin) * bsize + b - bmin;
		return data[n];
	}
	else {
		std::cout << "Error in index" << std::endl;
		return data[0];
	}
}

template<class T>
T& Matrix<T>::operator()(int a, int b, int c)
{
	if (a >= amin&&a <= amax&&b >= bmin&&b <= bmax&&c >= cmin&&c <= cmax) {
		int bsize = bmax - bmin + 1;
		int csize = cmax - cmin + 1;
		int n = (a - amin) * bsize * csize + (b - bmin) * csize + c - cmin;
		return data[n];
	}
	else
	{
		std::cout << "Error in index" << std::endl;
		return data[0];
	}
}

template<class T>
T& Matrix<T>::operator()(int a, int b, int c, int d) 
{
	if (a >= amin&&a <= amax&&b >= bmin&&b <= bmax&&c >= cmin&&c <= cmax&&d >= dmin&&d <= dmax) {
		int bsize = bmax - bmin + 1;
		int csize = cmax - cmin + 1;
		int dsize = dmax - dmin + 1;
		int n = (a - amin) * bsize * csize * dsize + (b - bmin) * csize * dsize
			+ (c - cmin) * dsize + d - dmin;
		return data[n];
	}
	else
	{
		std::cout << "Error in index" << std::endl;
		return data[0];
	}
}


template<class T>
const T& Matrix<T>::operator()(int a) const
{
	if (a >= amin&&a <= amax) {
		int n = a - amin;
		return data[n];
	}
	else {
		std::cout << "Error in index" << std::endl;
		return data[0];
	}

}

template<class T>
const T& Matrix<T>::operator()(int a, int b) const
{
	if (a >= amin&&a <= amax&&b >= bmin&&b <= bmax) {
		int bsize = bmax - bmin + 1;
		int n = (a - amin) * bsize + b - bmin;
		return data[n];
	}
	else {
		std::cout << "Error in index" << std::endl;
		return data[0];
	}
}

template<class T>
const T& Matrix<T>::operator()(int a, int b, int c) const
{
	if (a >= amin&&a <= amax&&b >= bmin&&b <= bmax&&c >= cmin&&c <= cmax) {
		int bsize = bmax - bmin + 1;
		int csize = cmax - cmin + 1;
		int n = (a - amin) * bsize * csize + (b - bmin) * csize + c - cmin;
		return data[n];
	}
	else
	{
		std::cout << "Error in index" << std::endl;
		return data[0];
	}
}

template<class T>
const T& Matrix<T>::operator()(int a, int b, int c, int d) const
{
	if (a >= amin&&a <= amax&&b >= bmin&&b <= bmax&&c >= cmin&&c <= cmax&&d >= dmin&&d <= dmax) {
		int bsize = bmax - bmin + 1;
		int csize = cmax - cmin + 1;
		int dsize = dmax - dmin + 1;
		int n = (a - amin) * bsize * csize * dsize + (b - bmin) * csize * dsize
			+ (c - cmin) * dsize + d - dmin;
		return data[n];
	}
	else
	{
		std::cout << "Error in index" << std::endl;
		return data[0];
	}
}