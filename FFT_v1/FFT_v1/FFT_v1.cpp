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
//ofstream ofs("Results.txt", ios::out);

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
	int p2=0, p3=0, p5=0;
	void ini_Array(); //initialize bitArray
	void BitReverse();
	void Butterfly();
	void Butterfly2();
	void Butterfly3();
	void Butterfly5();
public:
	void fft(int pow2, int pow3, int pow5);
	void fft2(int Num);
	void getX();
};

void FFT::fft2(int Num){
	while (Num % 2 == 0){
		p2++;
		Num /= 2;
	}
	while (Num % 3 == 0){
		p3++;
		Num /= 3;
	}
	while (Num % 5 == 0){
		p5++;
		Num /= 5;
	}
	BitReverse();
	Butterfly();
}

void FFT::fft(int pow2, int pow3, int pow5){
	p2 = pow2;
	p3 = pow3;
	p5 = pow5;
	BitReverse();
	Butterfly();
}

void FFT::Butterfly(){
	int k, p, q, r, s, t, c, m=1;
	Complex w_N, w_2N, w_3N, w_4N, tmp1, tmp2, tmp3, tmp4;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	double thetaN, theta;
	for (c = 0; c < p5; c++){
		m *= 5;
		thetaN = -2.0*M_PI / 5;
		w_N.Real = cos(thetaN);
		w_N.Imag = sin(thetaN);
		w_2N.Real = cos(2.0 * thetaN);
		w_2N.Imag = sin(2.0 * thetaN);
		w_3N.Real = cos(3.0 * thetaN);
		w_3N.Imag = sin(3.0 * thetaN);
		w_4N.Real = cos(4.0 * thetaN);
		w_4N.Imag = sin(4.0 * thetaN);
		for (k = 0; k < m / 5; k++)
		{
			theta = -2.0*k*M_PI / m;
			for (p = k; p<N; p += m)
			{
				q = p + m / 5;
				r = q + m / 5;
				s = r + m / 5;
				t = s + m / 5;
				// (Complex)tmp1 = (Complex)w * (Complex)X[q]
				tmp1.Real = cos(theta)*X[q].Real - sin(theta)*X[q].Imag;
				tmp1.Imag = cos(theta)*X[q].Imag + sin(theta)*X[q].Real;
				// (Complex)tmp2 = (Complex)w^2 * (Complex)X[r]
				tmp2.Real = cos(2 * theta)*X[r].Real - sin(2 * theta)*X[r].Imag;
				tmp2.Imag = cos(2 * theta)*X[r].Imag + sin(2 * theta)*X[r].Real;
				// (Complex)tmp3 = (Complex)w^3 * (Complex)X[s]
				tmp3.Real = cos(3 * theta)*X[s].Real - sin(3 * theta)*X[s].Imag;
				tmp3.Imag = cos(3 * theta)*X[s].Imag + sin(3 * theta)*X[s].Real;
				// (Complex)tmp4 = (Complex)w^4 * (Complex)X[t]
				tmp4.Real = cos(4 * theta)*X[t].Real - sin(4 * theta)*X[t].Imag;
				tmp4.Imag = cos(4 * theta)*X[t].Imag + sin(4 * theta)*X[t].Real;
				// (Complex)X[q] = (Complex)X[p] + (Complex)tmp1*w_N + (Complex)tmp2*w_2N + (Complex)tmp3*w_3N + (Complex)tmp4*w_4N
				X[q].Real = X[p].Real + (w_N.Real*tmp1.Real - w_N.Imag*tmp1.Imag) + (w_2N.Real*tmp2.Real - w_2N.Imag*tmp2.Imag) + (w_3N.Real*tmp3.Real - w_3N.Imag*tmp3.Imag) + (w_4N.Real*tmp4.Real - w_4N.Imag*tmp4.Imag);
				X[q].Imag = X[p].Imag + (w_N.Real*tmp1.Imag + w_N.Imag*tmp1.Real) + (w_2N.Real*tmp2.Imag + w_2N.Imag*tmp2.Real) + (w_3N.Real*tmp3.Imag + w_3N.Imag*tmp3.Real) + (w_4N.Real*tmp4.Imag + w_4N.Imag*tmp4.Real);
				// (Complex)X[r] = (Complex)X[p] + (Complex)tmp1*w_2N + (Complex)tmp2*w_4N + (Complex)tmp3*w_N + (Complex)tmp4*w_3N
				X[r].Real = X[p].Real + (w_2N.Real*tmp1.Real - w_2N.Imag*tmp1.Imag) + (w_4N.Real*tmp2.Real - w_4N.Imag*tmp2.Imag) + (w_N.Real*tmp3.Real - w_N.Imag*tmp3.Imag) + (w_3N.Real*tmp4.Real - w_3N.Imag*tmp4.Imag);
				X[r].Imag = X[p].Imag + (w_2N.Real*tmp1.Imag + w_2N.Imag*tmp1.Real) + (w_4N.Real*tmp2.Imag + w_4N.Imag*tmp2.Real) + (w_N.Real*tmp3.Imag + w_N.Imag*tmp3.Real) + (w_3N.Real*tmp4.Imag + w_3N.Imag*tmp4.Real);
				// (Complex)X[s] = (Complex)X[p] + (Complex)tmp1*w_3N + (Complex)tmp2*w_N + (Complex)tmp3*w_4N + (Complex)tmp4*w_2N
				X[s].Real = X[p].Real + (w_3N.Real*tmp1.Real - w_3N.Imag*tmp1.Imag) + (w_N.Real*tmp2.Real - w_N.Imag*tmp2.Imag) + (w_4N.Real*tmp3.Real - w_4N.Imag*tmp3.Imag) + (w_2N.Real*tmp4.Real - w_2N.Imag*tmp4.Imag);
				X[s].Imag = X[p].Imag + (w_3N.Real*tmp1.Imag + w_3N.Imag*tmp1.Real) + (w_N.Real*tmp2.Imag + w_N.Imag*tmp2.Real) + (w_4N.Real*tmp3.Imag + w_4N.Imag*tmp3.Real) + (w_2N.Real*tmp4.Imag + w_2N.Imag*tmp4.Real);
				// (Complex)X[t] = (Complex)X[p] + (Complex)tmp1*w_4N + (Complex)tmp2*w_3N + (Complex)tmp3*w_2N + (Complex)tmp4*w_N
				X[t].Real = X[p].Real + (w_4N.Real*tmp1.Real - w_4N.Imag*tmp1.Imag) + (w_3N.Real*tmp2.Real - w_3N.Imag*tmp2.Imag) + (w_2N.Real*tmp3.Real - w_2N.Imag*tmp3.Imag) + (w_N.Real*tmp4.Real - w_N.Imag*tmp4.Imag);
				X[t].Imag = X[p].Imag + (w_4N.Real*tmp1.Imag + w_4N.Imag*tmp1.Real) + (w_3N.Real*tmp2.Imag + w_3N.Imag*tmp2.Real) + (w_2N.Real*tmp3.Imag + w_2N.Imag*tmp3.Real) + (w_N.Real*tmp4.Imag + w_N.Imag*tmp4.Real);
				// (Complex)X[p] = (Complex)X[p] + (Complex)tmp1 + (Complex)tmp2
				X[p].Real += tmp1.Real + tmp2.Real + tmp3.Real + tmp4.Real;
				X[p].Imag += tmp1.Imag + tmp2.Imag + tmp3.Imag + tmp4.Imag;
			}
		}
	}
	for (c = 0; c < p3; c++){
		m *= 3;
		thetaN = -2.0*M_PI / 3;
		w_N.Real = cos(thetaN);
		w_N.Imag = sin(thetaN);
		w_2N.Real = cos(2 * thetaN);
		w_2N.Imag = sin(2 * thetaN);
		for (k = 0; k < m / 3; k++)
		{
			theta = -2.0*k*M_PI / m;
			for (p = k; p<N; p += m)
			{
				q = p + m / 3;
				r = q + m / 3;
				// (Complex)tmp1 = (Complex)w * (Complex)X[q]
				tmp1.Real = cos(theta)*X[q].Real - sin(theta)*X[q].Imag;
				tmp1.Imag = cos(theta)*X[q].Imag + sin(theta)*X[q].Real;
				// (Complex)tmp2 = (Complex)w^2 * (Complex)X[r]
				tmp2.Real = cos(2 * theta)*X[r].Real - sin(2 * theta)*X[r].Imag;
				tmp2.Imag = cos(2 * theta)*X[r].Imag + sin(2 * theta)*X[r].Real;
				// (Complex)X[r] = (Complex)X[p] + (Complex)tmp1*w_N^2 + (Complex)tmp2*w_N^4
				//							   = (Complex)X[p] + (Complex)tmp1*w_2N + (Complex)tmp2*w_N
				X[r].Real = X[p].Real + (w_2N.Real*tmp1.Real - w_2N.Imag*tmp1.Imag) + (w_N.Real*tmp2.Real - w_N.Imag*tmp2.Imag);
				X[r].Imag = X[p].Imag + (w_2N.Real*tmp1.Imag + w_2N.Imag*tmp1.Real) + (w_N.Real*tmp2.Imag + w_N.Imag*tmp2.Real);
				// (Complex)X[q] = (Complex)X[p] + (Complex)tmp1*w_N + (Complex)tmp2*w_N^2
				//                             = (Complex)X[p] + (Complex)tmp1*w_N + (Complex)tmp2*w_2N
				X[q].Real = X[p].Real + (w_N.Real*tmp1.Real - w_N.Imag*tmp1.Imag) + (w_2N.Real*tmp2.Real - w_2N.Imag*tmp2.Imag);
				X[q].Imag = X[p].Imag + (w_N.Real*tmp1.Imag + w_N.Imag*tmp1.Real) + (w_2N.Real*tmp2.Imag + w_2N.Imag*tmp2.Real);
				// (Complex)X[p] = (Complex)X[p] + (Complex)tmp1 + (Complex)tmp2
				X[p].Real += tmp1.Real + tmp2.Real;
				X[p].Imag += tmp1.Imag + tmp2.Imag;
			}
		}
	}
	for (c = 0; c < p2; c++){
		m *= 2;
		thetaN = -2.0*M_PI / 2;
		w_N.Real = cos(thetaN);
		w_N.Imag = sin(thetaN);
		for (k = 0; k < m / 2; k++)
		{
			theta = -2.0*k*M_PI / m;
			for (p = k; p<N; p += m)
			{
				q = p + m / 2;
				// (Complex)tmp = (Complex)w * (Complex)X[q]
				tmp1.Real = cos(theta)*X[q].Real - sin(theta)*X[q].Imag;
				tmp1.Imag = cos(theta)*X[q].Imag + sin(theta)*X[q].Real;
				// (Complex)X[q] = (Complex)X[p] - (Complex)tmp
				X[q].Real = X[p].Real - tmp1.Real;
				X[q].Imag = X[p].Imag - tmp1.Imag;
				// (Complex)X[p] = (Complex)X[p] + (Complex)tmp
				X[p].Real += tmp1.Real;
				X[p].Imag += tmp1.Imag;
			}
		}
	}
}

void FFT::Butterfly2(){
	int k, p, q, m = 2;
	Complex w, w_N, tmp;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	while (m <= N)
	{
		for (k = 0; k<m / 2; k++)
		{
			w.Real = cos(-2.0*k*M_PI / m);
			w.Imag = sin(-2.0*k*M_PI / m);
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
	int k, p, q, r, m = 3;
	Complex w, w_N, w_2N, tmp1, tmp2;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	w_N.Real = cos(-2.0*M_PI / 3);
	w_N.Imag = sin(-2.0*M_PI / 3);
	w_2N.Real = cos(-2.0*2*M_PI / 3);
	w_2N.Imag = sin(-2.0*2*M_PI / 3);
	while (m <= N)
	{
		for (k = 0; k<m / 3; k++)
		{
			w.Real = cos(-2.0*k*M_PI / m);
			w.Imag = sin(-2.0*k*M_PI / m);
			for (p = k; p<N; p += m)
			{
				q = p + m / 3;
				r = q + m / 3;
				// (Complex)tmp1 = (Complex)w * (Complex)X[q]
				tmp1.Real = w.Real*X[q].Real - w.Imag*X[q].Imag;
				tmp1.Imag = w.Real*X[q].Imag + w.Imag*X[q].Real;
				// (Complex)tmp2 = (Complex)w^2 * (Complex)X[r]
				tmp2.Real = (w.Real *w.Real - w.Imag * w.Imag)*X[r].Real - (2 * w.Real*w.Imag)*X[r].Imag;
				tmp2.Imag = (w.Real *w.Real - w.Imag * w.Imag)*X[r].Imag + (2 * w.Real*w.Imag)*X[r].Real;
				// (Complex)X[r] = (Complex)X[p] + (Complex)tmp1*w_N^2 + (Complex)tmp2*w_N^4
				//							   = (Complex)X[p] + (Complex)tmp1*w_2N + (Complex)tmp2*w_N
				X[r].Real = X[p].Real + (w_2N.Real*tmp1.Real - w_2N.Imag*tmp1.Imag) + (w_N.Real*tmp2.Real - w_N.Imag*tmp2.Imag);
				X[r].Imag = X[p].Imag + (w_2N.Real*tmp1.Imag + w_2N.Imag*tmp1.Real) + (w_N.Real*tmp2.Imag + w_N.Imag*tmp2.Real);
				// (Complex)X[q] = (Complex)X[p] + (Complex)tmp1*w_N + (Complex)tmp2*w_N^2
				//                             = (Complex)X[p] + (Complex)tmp1*w_N + (Complex)tmp2*w_2N
				X[q].Real = X[p].Real + (w_N.Real*tmp1.Real - w_N.Imag*tmp1.Imag) + (w_2N.Real*tmp2.Real - w_2N.Imag*tmp2.Imag);
				X[q].Imag = X[p].Imag + (w_N.Real*tmp1.Imag + w_N.Imag*tmp1.Real) + (w_2N.Real*tmp2.Imag + w_2N.Imag*tmp2.Real);
				// (Complex)X[p] = (Complex)X[p] + (Complex)tmp1 + (Complex)tmp2
				X[p].Real += tmp1.Real + tmp2.Real;
				X[p].Imag += tmp1.Imag + tmp2.Imag;
			}
		}
		m *=3;
	}
}

void FFT::Butterfly5(){
	int k, p, q, r, s, t, m = 5;
	Complex w, w_N, w_2N, w_3N, w_4N, tmp1, tmp2, tmp3, tmp4;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	double thetaN = -2.0*M_PI / 5,theta;
	w_N.Real = cos(thetaN);
	w_N.Imag = sin(thetaN);
	w_2N.Real = cos(2.0 * thetaN);
	w_2N.Imag = sin(2.0 * thetaN);
	w_3N.Real = cos(3.0 * thetaN);
	w_3N.Imag = sin(3.0 * thetaN);
	w_4N.Real = cos(4.0 * thetaN);
	w_4N.Imag = sin(4.0 * thetaN);

	while (m <= N)
	{
		for (k = 0; k<m / 5; k++)
		{
			theta = -2.0*k*M_PI / m;
			//w.Real = cos(2.0*k*M_PI / m);
			//w.Imag = sin(2.0*k*M_PI / m);
			for (p = k; p<N; p += m)
			{
				q = p + m / 5;
				r = q + m / 5;
				s = r + m / 5;
				t = s + m / 5;
				// (Complex)tmp1 = (Complex)w * (Complex)X[q]
				tmp1.Real = cos(theta)*X[q].Real - sin(theta)*X[q].Imag;
				tmp1.Imag = cos(theta)*X[q].Imag + sin(theta)*X[q].Real;
				// (Complex)tmp2 = (Complex)w^2 * (Complex)X[r]
				tmp2.Real = cos(2 * theta)*X[r].Real - sin(2 * theta)*X[r].Imag;
				tmp2.Imag = cos(2 * theta)*X[r].Imag + sin(2 * theta)*X[r].Real;
				// (Complex)tmp3 = (Complex)w^3 * (Complex)X[s]
				tmp3.Real = cos(3 * theta)*X[s].Real - sin(3 * theta)*X[s].Imag;
				tmp3.Imag = cos(3 * theta)*X[s].Imag + sin(3 * theta)*X[s].Real;
				// (Complex)tmp4 = (Complex)w^4 * (Complex)X[t]
				tmp4.Real = cos(4 * theta)*X[t].Real - sin(4 * theta)*X[t].Imag;
				tmp4.Imag = cos(4 * theta)*X[t].Imag + sin(4 * theta)*X[t].Real;
				// (Complex)X[q] = (Complex)X[p] + (Complex)tmp1*w_N + (Complex)tmp2*w_2N + (Complex)tmp3*w_3N + (Complex)tmp4*w_4N
				X[q].Real = X[p].Real + (w_N.Real*tmp1.Real - w_N.Imag*tmp1.Imag) + (w_2N.Real*tmp2.Real - w_2N.Imag*tmp2.Imag) + (w_3N.Real*tmp3.Real - w_3N.Imag*tmp3.Imag) + (w_4N.Real*tmp4.Real - w_4N.Imag*tmp4.Imag);
				X[q].Imag = X[p].Imag + (w_N.Real*tmp1.Imag + w_N.Imag*tmp1.Real) + (w_2N.Real*tmp2.Imag + w_2N.Imag*tmp2.Real) + (w_3N.Real*tmp3.Imag + w_3N.Imag*tmp3.Real) + (w_4N.Real*tmp4.Imag + w_4N.Imag*tmp4.Real);
				// (Complex)X[r] = (Complex)X[p] + (Complex)tmp1*w_2N + (Complex)tmp2*w_4N + (Complex)tmp3*w_N + (Complex)tmp4*w_3N
				X[r].Real = X[p].Real + (w_2N.Real*tmp1.Real - w_2N.Imag*tmp1.Imag) + (w_4N.Real*tmp2.Real - w_4N.Imag*tmp2.Imag) + (w_N.Real*tmp3.Real - w_N.Imag*tmp3.Imag) + (w_3N.Real*tmp4.Real - w_3N.Imag*tmp4.Imag);
				X[r].Imag = X[p].Imag + (w_2N.Real*tmp1.Imag + w_2N.Imag*tmp1.Real) + (w_4N.Real*tmp2.Imag + w_4N.Imag*tmp2.Real) + (w_N.Real*tmp3.Imag + w_N.Imag*tmp3.Real) + (w_3N.Real*tmp4.Imag + w_3N.Imag*tmp4.Real);
				// (Complex)X[s] = (Complex)X[p] + (Complex)tmp1*w_3N + (Complex)tmp2*w_N + (Complex)tmp3*w_4N + (Complex)tmp4*w_2N
				X[s].Real = X[p].Real + (w_3N.Real*tmp1.Real - w_3N.Imag*tmp1.Imag) + (w_N.Real*tmp2.Real - w_N.Imag*tmp2.Imag) + (w_4N.Real*tmp3.Real - w_4N.Imag*tmp3.Imag) + (w_2N.Real*tmp4.Real - w_2N.Imag*tmp4.Imag);
				X[s].Imag = X[p].Imag + (w_3N.Real*tmp1.Imag + w_3N.Imag*tmp1.Real) + (w_N.Real*tmp2.Imag + w_N.Imag*tmp2.Real) + (w_4N.Real*tmp3.Imag + w_4N.Imag*tmp3.Real) + (w_2N.Real*tmp4.Imag + w_2N.Imag*tmp4.Real);
				// (Complex)X[t] = (Complex)X[p] + (Complex)tmp1*w_4N + (Complex)tmp2*w_3N + (Complex)tmp3*w_2N + (Complex)tmp4*w_N
				X[t].Real = X[p].Real + (w_4N.Real*tmp1.Real - w_4N.Imag*tmp1.Imag) + (w_3N.Real*tmp2.Real - w_3N.Imag*tmp2.Imag) + (w_2N.Real*tmp3.Real - w_2N.Imag*tmp3.Imag) + (w_N.Real*tmp4.Real - w_N.Imag*tmp4.Imag);
				X[t].Imag = X[p].Imag + (w_4N.Real*tmp1.Imag + w_4N.Imag*tmp1.Real) + (w_3N.Real*tmp2.Imag + w_3N.Imag*tmp2.Real) + (w_2N.Real*tmp3.Imag + w_2N.Imag*tmp3.Real) + (w_N.Real*tmp4.Imag + w_N.Imag*tmp4.Real);
				// (Complex)X[p] = (Complex)X[p] + (Complex)tmp1 + (Complex)tmp2
				X[p].Real += tmp1.Real + tmp2.Real + tmp3.Real + tmp4.Real;
				X[p].Imag += tmp1.Imag + tmp2.Imag + tmp3.Imag + tmp4.Imag;
			}
		}
		m *= 5;
	}
}

void FFT::ini_Array(){
	int i = 0;
	int sizen = pow(2, p2)*pow(3, p3)*pow(5, p5);
	X = new Complex[sizen];
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

void FFT::getX(){
	int i;
	for (i = 0; i < pow(2, p2)*pow(3, p3)*pow(5, p5); i++){
		X[i].data();
		cout << endl;
	}
}

int main()
{
	clock_t t1, t2;
	FFT t;
	int a, b, c, N, Icase;
	cout << "Input (1) power or (2) number? ";
	cin >> Icase;
	switch (Icase){
	case 1:
		cout << "power of 2: ";
		cin >> a;
		cout << "power of 3: ";
		cin >> b;
		cout << "power of 5: ";
		cin >> c;
		t1 = clock();
		t.fft(a, b, c);
		t2 = clock();
		cout << "FFT of N = 0 ~ " << pow(2, a)*pow(3, b)*pow(5, c) << endl;
		break;
	case 2:
		cout << "Input the number: ";
		cin >> N;
		t1 = clock();
		t.fft2(N);
		t2 = clock();
		cout << "FFT of N = 0 ~ " << N << endl;
		break;
	default:
		break;
	}
	
	/*
	t1 = clock();
	t.fft(a, b, c);
	//t.fft2(360);
	t2 = clock();
	*/
	
	t.getX();
	printf("time = %f\n", (t2 - t1) / (double)(CLOCKS_PER_SEC));
	system("pause");
	return 0;
}

