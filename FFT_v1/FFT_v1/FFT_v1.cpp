// FFT_v1.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <Windows.h>
#include <iostream>
#include <cmath>
#include <fstream>
#define DEBUG 0
using namespace std;
const int MaxNum = 10000000;//2147483647;
const double M_PI = 3.14159265358979323846;
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
	int bitArray[100];
	int checkArray[100];
	int p2, p3, p5;
public:
	void ini_Array(); //initialize bitArray and checkArray
	void swap(Complex &a, Complex &b);
	void BitReverse();
	void Butterfly();
	void printX(int num); 
	void fft(int pow2, int pow3, int pow5);
};

void FFT::fft(int pow2, int pow3, int pow5){
	p2 = pow2;
	p3 = pow3;
	p5 = pow5;
	BitReverse();
	cout << endl;
	int n = pow(2, pow2)*pow(3, pow3)*pow(5, pow5);
	printX(n);
	cout << endl;
}

void FFT::Butterfly(){

}

void FFT::ini_Array(){
	int i;
	X = new Complex[MaxNum];
	for (i = 0; i < p2; i++){
		bitArray[i] = 1;
	}
	for (i = 0; i < p3; i++){
		bitArray[i + p2] = 2;
	}
	for (i = 0; i < p5; i++){
		bitArray[i + p2 + p3] = 4;
	}
	int n = pow(2, p2)*pow(3, p3)*pow(5, p5);
	X[n-1].Real = (double)(n-1)*1.0;
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
		printf("%d -> %d\n", p, q);
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
	if (DEBUG){
		cout << endl << "DEBUG" << endl;
		cout << "(P2, P3, P5) = " << p2 << ", " << p3 << ", " << p5 << endl;
		for (int i = 0; i < p2 + p3 + p5; i++){
			cout << bitArray[i] << " ";
		}
		cout << endl << "N = " << N << endl;
		cout << "m = " << m << endl;
		cout << "sum = " << sum << endl;
	}
}

void FFT::swap(Complex &a, Complex &b)
{
	Complex tmp = a;
	a = b;
	b = tmp;
}

void FFT::printX(int num){
	for (int i = 0; i < num; i++){
		X[i].data();
		cout << endl;
	}
}

int main()
{
	FFT t;
	t.fft(2,1,0);
	cin.get();
	return 0;
}

