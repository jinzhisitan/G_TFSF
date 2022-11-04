#pragma once

//PURPOSE:
//	 ��̬����1-4ά������ά�ȵ����飬����ÿһά�ȵ������޿���Ϊ������
//NOTES:
//	1.���������Ĭ�Ϲ��캯��Matrix()������Ҫ��һ��������Ӧά�ȵ�Allocate()����Ϊ�������ռ�͸�ֵ��
//  2.�����ڴ����з�ʽ���������ά���ݡ�
//	3.��Ҫ�����ع�����operator()���������е�ĳ��ֵ��
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
	int amin, amax; // ��1ά�ȵ����޺�����,��Ϊ������amax>=amin����
	int	bmin, bmax; // ��2ά�ȵ����޺�����,��Ϊ������bmax>=bmin����
	int	cmin, cmax; // ��3ά�ȵ����޺�����,��Ϊ������cmax>=cmin����
	int dmin, dmax; // ��4ά�ȵ����޺�����,��Ϊ������dmax>=dmin����

	T* data;   //��̬���������
public:
	//Ĭ�Ϲ��캯����������Ĭ�Ϲ��캯��������Ҫ��һ��������Ӧά�ȵ�Allocate()��������ռ�͸�ֵ��
	Matrix();  
	//Ϊ1ά���������ռ�͸�ֵ��
    //a1:���ޣ�a2:����; val: �����ֵ��
	Matrix(int a1, int a2, T val=0);
	Matrix(int a1, int a2, int b1, int b2, T val = 0);
	Matrix(int a1, int a2, int b1, int b2, int c1, int c2, T val = 0);
	//Ϊ4ά���������ռ�͸�ֵ��
    //a1:��1ά���ޣ�a2:��1ά����;  b1:��2ά���ޣ�b2:��2ά����;
	//c1:��3ά���ޣ�c2:��3ά����;  d1:��4ά���ޣ�d2:��4ά����;
	// val: �����ֵ��
	Matrix(int a1, int a2, int b1, int b2, int c1, int c2, int d1, int d2, T val = 0);
	//copy constructor
	Matrix(const Matrix &orig);
	//�ع�����ֵ�������
	Matrix & operator=(const Matrix &orig);
	//deconstructor
	~Matrix();

	//1ά�������ռ䡣
	//a1:���ޣ�a2:����; val: �����ֵ��
	void Allocate(int a1, int a2, T val=0);
	void Allocate(int a1, int a2, int b1, int b2, T val=0);
	void Allocate(int a1, int a2, int b1, int b2, int c1, int c2, T val=0);
	//Ϊ4ά���������ռ�͸�ֵ��
	//a1:��1ά���ޣ�a2:��1ά����;  b1:��2ά���ޣ�b2:��2ά����;
	//c1:��3ά���ޣ�c2:��3ά����;  d1:��4ά���ޣ�d2:��4ά����;
	//val: �����ֵ��
	void Allocate(int a1, int a2, int b1, int b2, int c1, int c2, int d1, int d2, T val=0);
	//initial array to @c initval
	//@param initval
	void initArray(T val = 0);

	//overload the parents��
	//����1ά�����е�data[a-amin]ֵ
	T & operator()(int a);
	T & operator()(int a, int b);
	T & operator()(int a, int b, int c);
	//����4ά�����е�data[a-amin][b - amin][c - amin][d - amin]
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