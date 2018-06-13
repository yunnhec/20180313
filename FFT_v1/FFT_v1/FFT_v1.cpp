// FFT_v1.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
//#include "fft_v1.h"
#include <time.h>
#include <omp.h>
using namespace std;

#include <Windows.h>
#include <iostream>
#include <cmath>
const double M_PI = 3.141592653589;

struct Complex{
	double Real = 0.0;
	double Imag = 0.0;
	void data(){
		printf("%.4f + %.4fi", Real, Imag);
	}
};

class FFT{
private:
	int p2 = 0, p3 = 0, p5 = 0;
	int *bitArray;
	Complex *X;
	void ini_ArrayX(); //initialize bitArray for all case
	void BitReverse();
	void BitReverse_dct();
	void BitReverse_dst();
	void Butterfly();
public:
	void fft(int pow2, int pow3, int pow5);
	void fft2(int Num);
	void dct2(int Num);
	void dct2_def(int Num);
	void dst1(int Num);
	void dst1_def(int Num);
	void getX();
	void getX_dct();
	void getX_dst();
};

void FFT::dst1(int Num){
	Num = 2 * Num + 2;
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
	cout << "Num = " << Num << endl;
	//if (Num % 2 == 0 && Num % 3 == 0 && Num % 5 == 0)
	if (Num == 1)
	{
		BitReverse_dst();
		Butterfly();
	}
	else{
		cout << "The input cannot be solved by this program." << endl;
	}
}

void FFT::dct2(int Num){
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
	p2 += 2;
	BitReverse_dct();
	Butterfly();
}

void FFT::ini_ArrayX(){
	int i = 0;
	int sizen = pow(2, p2)*pow(3, p3)*pow(5, p5);
	bitArray = new int[p2 + p3 + p5];
	X = new Complex[sizen];
	for (i = 0; i < p2; i++)
		bitArray[i] = 1;
	for (i = 0; i < p3; i++)
		bitArray[i + p2] = 2;
	for (i = 0; i < p5; i++)
		bitArray[i + p2 + p3] = 4;
}

void FFT::BitReverse_dct(){
	ini_ArrayX();
	int m, p, q, k, c = 1;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	int sum = p2 + p3 + p5;
	m = N / (bitArray[sum - c] + 1);
	q = m;
	for (p = 1; p<N - 1; ++p)
	{
		if (q>N / 2 || q % 2 == 0){
			X[p].Real = 0.0;
		}
		else{
			X[p].Real = (q - 1) / 2.0;
		}
		k = m;
		while (q >= bitArray[sum - c] * k & k>0) {
			q = q - bitArray[sum - c] * k;
			k = k / (bitArray[sum - c - 1] + 1);
			c++;
		}
		c = 1;
		q = q + k;
	}
	X[N - 1].Real = 0.0;
}

void FFT::BitReverse_dst(){
	ini_ArrayX();
	int m, p, q, k, c = 1;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	int sum = p2 + p3 + p5;
	m = N / (bitArray[sum - c] + 1);
	q = m;
	for (p = 1; p<N - 1; ++p)
	{
		if (q == N / 2){
			X[p].Real = 0.0;
		}
		else if (q <= N / 2) {	//q is not zero
			X[p].Real = (double)q-1;
		}
		else if (q>N / 2) {
			X[p].Real = (double)(q+1-N);
		}
		k = m;
		while (q >= bitArray[sum - c] * k & k>0) {
			q = q - bitArray[sum - c] * k;
			k = k / (bitArray[sum - c - 1] + 1);
			c++;
		}
		c = 1;
		q = q + k;
	}
}

void FFT::dct2_def(int num){
	int n, k;
	double *A = new double[num];
	int *x = new int[num];
	//initialize
	for (k = 0; k < num; k++){
		A[k] = 0;
		x[k] = k;
	}
	for (k = 0; k < num; k++){
		for (n = 0; n < num; n++){
			A[k] = A[k] + x[n] * cos(M_PI / (2 * num)*k*(2 * n + 1));
		}
	}
	cout << "Calculate by definition: " << endl;
	for (k = 0; k < num; k++){
		printf("%.4f\n", A[k]);
	}
}

void FFT::dst1_def(int num){
	int n, k;
	double *A = new double[num];
	int *x = new int[num];
	//initialize
	for (k = 0; k < num; k++){
		A[k] = 0;
		x[k] = k;
	}
	for (k = 0; k < num; k++){
		for (n = 0; n < num; n++){
			A[k] = A[k] + x[n] * sin(M_PI / (num+1)*(k+1)*(n + 1));
		}
	}
	cout << "Calculate by definition: " << endl;
	for (k = 0; k < num; k++){
		printf("%.4f\n", A[k]);
	}
}

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
	int k, p, q, r, s, t, c, m = 1;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	double thetaN, theta;
	Complex w_N, w_2N, w_3N, w_4N, tmp1, tmp2, tmp3, tmp4;
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
			//#pragma omp for num_threads(4)
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
			//#pragma omp for num_threads(4)
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
			//#pragma omp for num_threads(4)
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

void FFT::BitReverse(){
	ini_ArrayX();
	int m, p, q, k, c = 1;
	int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	int sum = p2 + p3 + p5;
	m = N / (bitArray[sum - c] + 1);
	q = m;
	//#pragma omp for num_threads(4)
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
	cout << "FFT of N = 0 ~ " << pow(2, p2)*pow(3, p3)*pow(5, p5) << endl;

	for (i = 0; i < pow(2, p2)*pow(3, p3)*pow(5, p5); i++){
		X[i].data();
		cout << endl;
	}
	cout << endl;
}

void FFT::getX_dct(){
	int i;
	cout << "DCT-II of N = 0 ~ " << pow(2, p2-2)*pow(3, p3)*pow(5, p5)-1 << endl;

	for (i = 0; i < pow(2, p2-2)*pow(3, p3)*pow(5, p5); i++){
		printf("%.4f\n", X[i].Real);
	}
	cout << endl;
}

void FFT::getX_dst(){
	int i;
	cout << "DST-I of N = 0 ~ " << pow(2, p2)*pow(3, p3)*pow(5, p5)/2 - 1 << endl;

	for (i = 1; i < pow(2, p2)*pow(3, p3)*pow(5, p5)/2; i++){
		printf("%.4f\n", - X[i].Imag / 2.0);
	}
	cout << endl;
}

int main()
{
	clock_t t1, t2;
	FFT t;
	int input;
	cout << "N = ";
	cin >> input;
	t1 = clock();
	t.dst1(input);
	//t.dct2(input);
	t2 = clock();
	printf("time = %f\n", (t2 - t1) / (double)(CLOCKS_PER_SEC));
	system("pause");
	t.getX_dst();
	//t.getX_dct();
	system("pause");
	t.dst1_def(input);
	//t.dct2_def(input);
	system("pause");
	return 0;
}

