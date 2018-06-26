// FFT_Zp_v1.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <Windows.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#define DEBUG 10
using namespace std;

class FFT_v2{
private:
	int i, L1, L2, digit, Iter_Max = 1000;
	int P, Wn=1, invWn=1;
	int p2 = 0, p3 = 0, p5 = 0;
	int *X, *Y, *Z, *bitArray;
	void find_wn();
	void BitReverse2(string, string);
	void ButterflyZp();
	void iFFT();
	void GE();
	bool isPrime(int); //return true if input number is a prime
	void fact(int); //factorization of digit
public:
	void Multi(string, string);
	void print();
};

void FFT_v2::fact(int Num){
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
}

bool FFT_v2::isPrime(int number){
	for (i = 2; i*i <= number; i++)
	{
		if (number % i == 0)
			return false;
	}
	return true;
}

void FFT_v2::Multi(string in1, string in2){
	BitReverse2(in1, in2);
	find_wn();
}

void FFT_v2::find_wn(){
	cout << "digit = " << digit << endl;
	for (int k = 0; k < Iter_Max; k++){
		P = 1 + k*digit;
		if (isPrime(P) && P>max(L1, L2) * 81){
			// find a W such that W^digit % P = 1
			for (Wn = 2; Wn < P; Wn++){
				int tmp = Wn;
				//operate Wn^digit
				for (i = 2; i < P; i++){
					tmp = tmp*Wn % P;
					if (tmp == 1)
						break;
				}
				if (i==digit)
					break;
			}
			if (Wn < P)
				break;
		}
	}
	//Wn^-1 = Wn^(digit-1) % P
	for (i = 0; i<digit - 1; i++)
	{
		invWn = invWn*Wn %P;
	}
	if (DEBUG)
		cout << "(P, Wn, invWn) =" << P << " " << Wn<< " " << invWn << endl;
}

void FFT_v2::BitReverse2(string n1,string n2){
	//ini_Array();
	L1 = n1.length();
	L2 = n2.length();
	digit = L1 + L2;
	if (digit % 2 != 0){
		if (digit % 3 != 0 && digit % 5 != 0) //if 2,3,5 are not factor of digit, then digit++
			digit++;
	}
	X = new int[digit];
	Y = new int[digit];
	Z = new int[digit];

	//factorization of digit
	fact(digit);
	int sizen = pow(2, p2)*pow(3, p3)*pow(5, p5);
	bitArray = new int[p2 + p3 + p5];
	for (i = 0; i < p2; i++)
		bitArray[i] = 1;
	for (i = 0; i < p3; i++)
		bitArray[i + p2] = 2;
	for (i = 0; i < p5; i++)
		bitArray[i + p2 + p3] = 4;

	
	//BitReverse of 2,3,5
	int m, p, q, k, c = 1;
	// int N = pow(2, p2)*pow(3, p3)*pow(5, p5);
	int sum = p2 + p3 + p5;
	m = digit / (bitArray[sum - c] + 1);
	q = m;
	X[0] = n1[L1-1] - 48;
	Y[0] = n2[L2-1] - 48;
	for (p = 1; p<digit - 1; ++p)
	{
		//printf("%d -> %d\n", p, q);
		if (q>=L1){
			X[p] = 0;
		}
		else{
			X[p] = n1[L1 - 1 - q] - 48;
		}
		if (q >= L2){
			Y[p] = 0;
		}
		else{
			Y[p] = n2[L2 - 1 - q] - 48;
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
	X[digit - 1] = 0;
	Y[digit - 1] = 0;
	if (DEBUG){
		cout << "X = ";
		for (i = 0; i < digit; i++){
			cout << X[i] << " ";
		}
		cout << endl << "Y = ";
		for (i = 0; i < digit; i++){
			cout << Y[i] << " ";
		}
		cout << endl;
	}
}

void FFT_v2::GE(){
	int i, j, s, a, b;
	a = 1;
	for (i = 0; i<digit; ++i)
	{
		s = 0;
		b = 1;
		for (j = 0; j<digit; ++j)
		{
			s = (s + b*Y[j]) % P;
			b = (b*a) % P;
		}
		Z[i] = s;
		a = (a*Wn) % P;
	}
}
/*
void FFT_v2::ButterflyZp(){
	int k, p, q, r, s, t, c, m = 1,thet;
	int Ww, w_N, w_2N, w_3N, w_4N, tmp1, tmp2, tmp3, tmp4;
	/*for (c = 0; c < p5; c++){
		m *= 5;
		w_N = invWn;
		w_2N = invWn*invWn % P;
		w_3N = w_2N*invWn % P;
		w_4N = w_3N*invWn % P;
		for (k = 0; k < m / 5; k++)
		{
			//thetaN = -2.0*M_PI / 5;
			theta = -2.0*k*M_PI / m;
			//#pragma omp for num_threads(4)
			for (p = k; p<N; p += m)
			{
				q = p + m / 5;
				r = q + m / 5;
				s = r + m / 5;
				t = s + m / 5;
				// (Complex)tmp1 = (Complex)w * (Complex)X[q]
				
				// (Complex)tmp2 = (Complex)w^2 * (Complex)X[r]
				
				// (Complex)tmp3 = (Complex)w^3 * (Complex)X[s]

				// (Complex)tmp4 = (Complex)w^4 * (Complex)X[t]

				// (Complex)X[q] = (Complex)X[p] + (Complex)tmp1*w_N + (Complex)tmp2*w_2N + (Complex)tmp3*w_3N + (Complex)tmp4*w_4N

				// (Complex)X[r] = (Complex)X[p] + (Complex)tmp1*w_2N + (Complex)tmp2*w_4N + (Complex)tmp3*w_N + (Complex)tmp4*w_3N

				// (Complex)X[s] = (Complex)X[p] + (Complex)tmp1*w_3N + (Complex)tmp2*w_N + (Complex)tmp3*w_4N + (Complex)tmp4*w_2N
				
				// (Complex)X[t] = (Complex)X[p] + (Complex)tmp1*w_4N + (Complex)tmp2*w_3N + (Complex)tmp3*w_2N + (Complex)tmp4*w_N
				
				// (Complex)X[p] = (Complex)X[p] + (Complex)tmp1 + (Complex)tmp2
				
			}
		}
	}
	for (c = 0; c < p3; c++){
		m *= 3;
		thetaN = -2.0*M_PI / 3;
		w_N = cos(thetaN);
		w_N = sin(thetaN);
		w_2N = cos(2 * thetaN);
		w_2N = sin(2 * thetaN);
		for (k = 0; k < m / 3; k++)
		{
			theta = -2.0*k*M_PI / m;
			//#pragma omp for num_threads(4)
			for (p = k; p<N; p += m)
			{
				q = p + m / 3;
				r = q + m / 3;
				// (Complex)tmp1 = (Complex)w * (Complex)X[q]
				tmp1 = cos(theta)*X[q].Real - sin(theta)*X[q].Imag;
				tmp1 = cos(theta)*X[q].Imag + sin(theta)*X[q].Real;
				// (Complex)tmp2 = (Complex)w^2 * (Complex)X[r]
				tmp2 = cos(2 * theta)*X[r].Real - sin(2 * theta)*X[r].Imag;
				tmp2 = cos(2 * theta)*X[r].Imag + sin(2 * theta)*X[r].Real;
				// (Complex)X[r] = (Complex)X[p] + (Complex)tmp1*w_N^2 + (Complex)tmp2*w_N^4
				//							   = (Complex)X[p] + (Complex)tmp1*w_2N + (Complex)tmp2*w_N
				X[r] = X[p].Real + (w_2N.Real*tmp1.Real - w_2N.Imag*tmp1.Imag) + (w_N.Real*tmp2.Real - w_N.Imag*tmp2.Imag);
				X[r] = X[p].Imag + (w_2N.Real*tmp1.Imag + w_2N.Imag*tmp1.Real) + (w_N.Real*tmp2.Imag + w_N.Imag*tmp2.Real);
				// (Complex)X[q] = (Complex)X[p] + (Complex)tmp1*w_N + (Complex)tmp2*w_N^2
				//                             = (Complex)X[p] + (Complex)tmp1*w_N + (Complex)tmp2*w_2N
				X[q] = X[p].Real + (w_N.Real*tmp1.Real - w_N.Imag*tmp1.Imag) + (w_2N.Real*tmp2.Real - w_2N.Imag*tmp2.Imag);
				X[q] = X[p].Imag + (w_N.Real*tmp1.Imag + w_N.Imag*tmp1.Real) + (w_2N.Real*tmp2.Imag + w_2N.Imag*tmp2.Real);
				// (Complex)X[p] = (Complex)X[p] + (Complex)tmp1 + (Complex)tmp2
				X[p].Real += tmp1.Real + tmp2.Real;
				X[p].Imag += tmp1.Imag + tmp2.Imag;
			}
		}
	}

	
	for (c = 0; c < p2; c++){
		m *= 2;
		Ww = 1;
		for (k = 0; k < m / 2; k++)
		{
			for (p = k; p<digit; p += m)
			{ 
				q = p + m / 2;
				//tmp = X[q]*Ww
				tmp1 = pow(invWn, m);
				tmp1 = X[q] * tmp1 %P;
				//X[q] = X[p]+tmp*Wn
				X[q] = (X[p] + tmp1*invWn) % P;
				//X[p] = X[p]+tmp
				X[p] = (X[p] + tmp1) % P;
			}
			
		}
	}
}
*/
void FFT_v2::print(){
	cout << "X = ";
	for (i = 0; i < digit; i++){
		cout << X[i] << " ";
	}
	cout << endl << "Y = ";
	for (i = 0; i < digit; i++){
		cout << Y[i] << " ";
	}
	cout << endl << "Z = ";
	for (i = 0; i < digit; i++){
		cout << Z[i] << " ";
	}
	cout << endl;
}


int main()
{
	FFT_v2 test;
	string input1, input2;
	/*cout << "input 1= ";
	cin >> input1;
	cout << "input 2= ";
	cin >> input2;*/
	input1 = "1324"; 
	input2 = "1324";
	test.Multi(input1, input2);
	system("pause");
	test.print();
	system("pause");
	return 0;
}

