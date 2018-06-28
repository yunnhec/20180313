// FFT_Zp_v1.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <Windows.h>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#define DEBUG 0
using namespace std;

class FFT_v2{
private:
	int i, L1, L2, digit, Iter_Max = 1000;
	int result=0, P, Wn = 1, invWn = 1, p2 = 0, p3 = 0, p5 = 0;
	int *X, *Y, *T,*Z, *bitArray;
	void find_wn();
	void BitReverse2(string, string);
	void fftZp();
	void multi();
	void ifftZp();
	void conv(); //convert string elements to numbers
	bool isPrime(int); //return true if input number is a prime
	void fact(int); //factorization of digit
	int Mod(int, int);
	int Power(int, int,int);
public:
	void Multi(string, string);
	void print();
	void mulans();
};

void FFT_v2::Multi(string in1, string in2){
	BitReverse2(in1, in2);
	find_wn();
	fftZp();
	multi();
	ifftZp();
	conv();
}

void FFT_v2::conv(){
	T[0] = Z[0];
	for (i = 1; i < digit; i++){
		T[i] = Z[digit - i];
	}
}

void FFT_v2::fftZp(){
	int Wtmp = invWn; //determine fft or ifft
	int k, p, q, c, m = 2;
	int wk,wt, tmp1;
	for (c = 0; c < p2; c++){
		wt = Mod(Power(Wtmp, digit / m, P), P);
		for (k = 0; k < m / 2; k++){
			for (p = k; p < digit; p += m){
				q = p + m / 2;

				tmp1 = X[q];
				X[q] = X[p] + Power(wt, k + m / 2,P)*tmp1;
				X[q] = X[q] % P;
				X[p] = X[p] + Power(wt, k,P)*tmp1;
				X[p] = X[p] % P;

				tmp1 = Y[q];
				Y[q] = Y[p] + Power(wt, k + m / 2, P)*tmp1;
				Y[q] = Y[q] % P;
				Y[p] = Y[p] + Power(wt, k, P)*tmp1;
				Y[p] = Y[p] % P;
			}
		}
		m *= 2;
	}
}

void FFT_v2::ifftZp(){
	int Wtmp = Wn; //determine fft or ifft
	int k, p, q, c, m = 2;
	int wk, wt, tmp1;
	for (c = 0; c < p2; c++){
		wt = Mod(Power(Wtmp, digit / m, P), P);
		for (k = 0; k < m / 2; k++){
			for (p = k; p < digit; p += m){
				q = p + m / 2;
				tmp1 = Z[q];
				Z[q] = Z[p] + Power(wt, k + m / 2, P)*tmp1;
				Z[q] = Z[q] % P;
				Z[p] = Z[p] + Power(wt, k, P)*tmp1;
				Z[p] = Z[p] % P;
			}
		}
		m *= 2;
	}
	for (i = 0; i < digit; i++){
		while (Z[i] % digit != 0){
			Z[i] += P;
		}
			Z[i] = Mod(Z[i] / digit, P);
	}
}

void FFT_v2::multi(){
	T[0] = Mod(X[0] * Y[0], P);
	for (i = 1; i < digit; i++){
		T[i] = Mod(X[digit - i] * Y[digit - i], P);
	}
	//bit reverse of z (z to Z)
	int m, p, q, k, c = 1;
	int sum = p2 + p3 + p5;
	m = digit / (bitArray[sum - c] + 1);
	q = m;
	Z[0] = T[0];
	for (p = 1; p<digit - 1; ++p)
	{
		//printf("%d -> %d\n", p, q);
		Z[q] = T[p];
		k = m;
		while (q >= bitArray[sum - c] * k & k>0) {
			q = q - bitArray[sum - c] * k;
			k = k / (bitArray[sum - c - 1] + 1);
			c++;
		}
		c = 1;
		q = q + k;
	}
	Z[digit - 1] = T[digit - 1];
}

void FFT_v2::find_wn(){
	//cout << "digit = " << digit << endl;
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
	T = new int[digit];
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
	int sum = p2 + p3 + p5;
	m = digit / (bitArray[sum - c] + 1);
	q = m;
	X[0] = n1[L1-1] - 48;
	Y[0] = n2[L2-1] - 48;
	for (p = 1; p<digit - 1; ++p)
	{
		//printf("%d -> %d\n", p, q);
		if (q >= L1){
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
#if DEBUG
		cout << "X = ";
		for (i = 0; i < digit; i++) {
			cout << X[i] << " ";
		}
		cout << endl << "Y = ";
		for (i = 0; i < digit; i++) {
			cout << Y[i] << " ";
		}
		cout << endl;
#endif
}

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

int FFT_v2::Mod(int x, int prime) {
	int t = x%prime;
	if (t < 0)
		t += prime;
	return t;
}

int FFT_v2::Power(int num, int power,int prime){
	int t = 1;
	if (power == 0)
		return 1;
	else{
		t = num;
		for (int j = 1; j < power; j++){
			t = t*num;
			t = t%prime;
		}
		return t;
	}
}

void FFT_v2::print(){
	/*cout << "X = ";
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
	}*/
	cout << "iFFT(X*Y) = ";
	for (i = 0; i < digit; i++){
		cout << T[i] << " ";
	}
	cout << endl;
	
}

void FFT_v2::mulans(){
	for (i = 1; i < digit; i++){
		T[i] = Z[digit - i];
	}
	for (i = 0; i<digit - 1; i++)
	{
		T[i + 1] += T[i] / 10;
		T[i] = T[i] % 10;
	}
	for (i = digit - 1; i >= 0; i--)
	{
		cout << T[i];
	}
	cout << endl;
}

int main()
{
	FFT_v2 test;
	string input1, input2;
	cout << "digit of number must be radix-2" << endl;
	cout << "X = ";
	cin >> input1;
	cout << "Y = ";
	cin >> input2;
	test.Multi(input1, input2);
	test.print();
	system("pause");
	test.mulans();
	system("pause");
	return 0;
}

