// FFT_v1.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <Windows.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#define DEBUG 0
using namespace std;
const int MaxNum = 100000000;//2147483647;
const double M_PI = 3.141592653589;
//ofstream ofs("Output.txt", ios::out);

struct Complex{
	double Real = 0.0;
	double Imag = 0.0;
	void data(){
		printf("%.4f + %.4fi", Real, Imag);
		//ofs << Real << " + " << Imag << "i";
	}
};

class FFT{
private:
	Complex *X;
	int bitArray[50];
	int p2, p3, p5;
public:
	void ini_Array(); //initialize bitArray
	void swap(Complex &a, Complex &b);
	void BitReverse();
	void Butterfly2();
	void Butterfly3();
	void getX(); 
	void fft(int pow2, int pow3, int pow5);
};

void FFT::fft(int pow2, int pow3, int pow5){
	p2 = pow2;
	p3 = pow3;
	p5 = pow5;
	BitReverse();
}

void FFT::Butterfly2(){
	int k, p, q, m = 1;
	Complex w, w_N, tmp;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	while (m <= N)
	{
		for (k = 0; k<m / 2; k++)
		{
			w.Real = cos(2.0*k*M_PI / m); // -????????? 
			w.Imag = sin(2.0*k*M_PI / m);
			for (p = k; p<N; p += m)
			{
				q = p + m / 2;
				// (Complex)tmp = (Complex)w * (Complex)X[q]
				tmp.Real = w.Real*X[q].Real - w.Imag*X[q].Imag;
				tmp.Imag = w.Real*X[q].Imag + w.Imag*X[q].Real;
				// (Complex)X[q] = (Complex)X[p] - (Complex)tmp
				X[q].Real = X[p].Real - tmp.Real;
				X[q].Imag = X[p].Imag - tmp.Imag;
				// (Complex)X[p] = (Complex)X[p] + (Complex)tmp
				X[p].Real += tmp.Real;
				X[p].Imag += tmp.Imag;
			}
		}
		m <<= 1; // m*=2;
	}
}

void FFT::Butterfly3(){
	int k, p, q, r, m = 1;
	Complex w, w_N, tmp1, tmp2;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	w_N.Real = cos(2.0*M_PI / 5);
	w_N.Imag = sin(2.0*M_PI / 5);
	while (m <= N)
	{
		for (k = 0; k<m / 3; k++)
		{
			w.Real = cos(2.0*k*M_PI / m);
			w.Imag = sin(2.0*k*M_PI / m);
			for (p = k; p<N; p += m)
			{
				q = p + m / 3;
				r = q + m / 3;
				// (Complex)tmp1 = (Complex)w * (Complex)X[q]
				// (Complex)tmp2 = (Complex)w^2 * (Complex)X[r]

				// (Complex)X[r] = (Complex)X[p] +...

				// (Complex)X[q] = (Complex)X[p] +...

				// (Complex)X[p] = (Complex)X[p] + (Complex)tmp1 + (Complex)tmp2
				
			}
		}
		m *=3;
	}
}

void FFT::ini_Array(){
	int i = 0;
	X = new Complex[MaxNum];
	for (i = 0; i < p2; i++)
		bitArray[i] = 1;
	for (i = 0; i < p3; i++)
		bitArray[i + p2] = 2;
	for (i = 0; i < p5; i++)
		bitArray[i + p2 + p3] = 4;
}

void FFT::BitReverse(){
	ini_Array();
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	int sum = p2 + p3 + p5;
	int m, p, q, k, c = 1;
	m = N / (bitArray[sum - c] + 1);
	q = m;
	for (p = 1; p<N - 1; ++p)
	{
		//printf("%d -> %d\n", p, q);
		X[p].Real = q*1.0;
		k = m;
		while (q >= bitArray[sum - c] * k & k>0) {
			q = q - bitArray[sum - c] * k;
			k = k / (bitArray[sum - c - 1] + 1);
			c++;
		}
		c = 1;
		q = q + k;
	}
	X[N - 1].Real = (double)(N - 1)*1.0;
}

void FFT::swap(Complex &a, Complex &b)
{
	Complex tmp = a;
	a = b;
	b = tmp;
}

void FFT::getX(){
	for (int i = 0; i < pow(2, p2)*pow(3, p3)*pow(5, p5); i++){
		X[i].data();
		cout << endl;
	}
}

int main()
{
	clock_t t1, t2;
	FFT t;
	cout << "BitReverse" << endl;
	t1 = clock();
	t.fft(4,0,0);
	t2 = clock();
	t.getX();
	t.Butterfly2();
	cout << "Butterfly" << endl;
	t.getX();
	printf("time = %f\n", (t2 - t1) / (double)(CLOCKS_PER_SEC));

	system("pause");
	return 0;
}

