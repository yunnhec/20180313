// FFT_v1.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <Windows.h>
#include <iostream>
#include <math.h>
#include <fstream>
#define DEBUG 0
using namespace std;

struct Complex{
	double Real = 0.0;
	double Imag = 0.0;
	void data(){
		printf("%d + i%d", Real, Imag);
	}
};

class FFT{
private:
	Complex A;
	int bitArray[100];
	int checkArray[100];
public:
	int BitReverse2(int num2);
	int BitReverse3(int num3);
	int BitReverse5(int num5);
	int BitReverse(int pow2, int pow3, int pow5);
	void ini_Array(int p2, int p3, int p5); //initialize bitArray and checkArray
	void swap(Complex &a, Complex &b);
	void test();
};

void FFT::test(){
	//BitReverse(1, 1, 1);
	Complex a;
	a.Imag = 1.5;
	a.Real = 1.2;
	a.data();
}

void FFT::ini_Array(int p2,int p3,int p5){
	int i;
	for (i = 0; i < p2; i++){
		bitArray[i] = 1;
		checkArray[i] = 0;
	}
	for (i = 0; i < p3; i++){
		bitArray[i + p2] = 2;
		checkArray[i] = 0;
	}
	for (i = 0; i < p5; i++){
		bitArray[i + p2 + p3] = 4;
		checkArray[i] = 0;
	}
}

int FFT::BitReverse(int pow2, int pow3, int pow5){
	ini_Array(pow2, pow3, pow5);
	int N = pow(2, pow2)*pow(3, pow3)*pow(5, pow5);
	int sum = pow2 + pow3 + pow5;
	int m, p, q, k, c=1;
	m = N / (bitArray[sum - c] + 1);
	q = m;
	for (p = 1; p<N - 1; ++p)
	{
		printf("%d -> %d\n", p, q);
		if (checkArray[p] == 0 && checkArray[q]==0)
		{
			//swap p and q
			checkArray[p] = 1;
			checkArray[q] = 1;
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
	if (DEBUG){
		cout << endl << "DEBUG" << endl;
		cout << "(P2, P3, P5) = " << pow2 << ", " << pow3 << ", " << pow5 << endl;
		for (int i = 0; i < pow2 + pow3 + pow5; i++){
			cout << bitArray[i] << " ";
		}
		cout << endl << "N = " << N << endl;
		cout << "m = " << m << endl;
		cout << "sum = " << sum << endl;
	}
	return 0;
}

int FFT::BitReverse2(int num2){
	int m, p, q, k;
	m = num2 / 2;                        // Bit-Reverse 每次要進位的數字 
	q = m;							// p = 1, q = m (第一個要交換的) 
	for (p = 1; p<num2 - 1; ++p)
	{
		printf("%d <-> %d\n", p, q);
		if (p < q)
		{
			//swap p and q
		}
		k = m;						// k, 用來檢查第 log_2 k + 1 位是不是1 
		while (q >= k & k > 0)		// q >=k 第 (log_2 k + 1)位是1,  
		{
			q = q - k;				// 1->0
			k = k / 2;				// 檢查下一位 
		}
		q = q + k;
	}
	return 0;
}

int FFT::BitReverse3(int num3){
	int m, p, q, k;
	m = num3 / 3;
	q = m;
	for (p = 1; p < num3-1; p++){
		printf("%d <-> %d\n", p, q);
		if (p < q)
		{
			//swap p and q
		}
		k = m;
		while (q >= 2*k & k > 0)		// q >=k 第 (log_3 k + 1)位是2,  
		{
			q = q - 2*k;				// 2->0
			k = k / 3;				// 檢查下一位 
		}
		q = q + k;
	}
	return 0;
}

int FFT::BitReverse5(int num5){
	int m, p, q, k;
	m = num5 / 5;
	q = m;
	for (p = 1; p < num5 - 1; p++){
		printf("%d <-> %d\n", p, q);
		if (p < q)
		{
			//swap p and q
		}
		k = m;
		while (q >= 4 * k & k > 0)		// q >=k 第 (log_5 k + 1)位是4,  
		{
			q = q - 4 * k;				// 4->0
			k = k / 5;				// 檢查下一位 
		}
		q = q + k;
	}
	return 0;
}

void FFT::swap(Complex &a, Complex &b)
{
	Complex tmp;
	tmp.Real = a.Real;
	tmp.Imag = a.Imag;
	a.Real = b.Real;
	a.Imag = b.Imag;
	b.Real = tmp.Real;
	b.Imag = tmp.Imag;
}

int main()
{
	FFT t;
	t.test();
	cin.get();
	return 0;
}

