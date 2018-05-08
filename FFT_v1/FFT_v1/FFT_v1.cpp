// FFT_v1.cpp : �w�q�D���x���ε{�����i�J�I�C
//

#include "stdafx.h"
#include <Windows.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include "Complex_1.h"
using namespace std;

class FFT{
private:
	Complex<double> *in, *out;
public:
	int BitReverse2(int num2);
	int BitReverse3(int num3);
	int BitReverse5(int num5);
	void swap(int &a, int &b);
	void test();
};

void FFT::test(){
	int x = 128, y = 1;
	cout << "(x,y) = " << x << ", " << y << endl;
	swap(x, y);
	cout << "(x,y) = " << x << ", " << y << endl;
}

int FFT::BitReverse2(int num2){
	int m, p, q, k;
	m = num2 / 2;                        // Bit-Reverse �C���n�i�쪺�Ʀr 
	q = m;							// p = 1, q = m (�Ĥ@�ӭn�洫��) 
	for (p = 1; p<num2 - 1; ++p)
	{
		printf("%d <-> %d\n", p, q);
		if (p < q)
		{
			//swap p and q
		}
		k = m;						// k, �Ψ��ˬd�� log_2 k + 1 ��O���O1 
		while (q >= k & k > 0)		// q >=k �� (log_2 k + 1)��O1,  
		{
			q = q - k;				// 1->0
			k = k / 2;				// �ˬd�U�@�� 
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
		while (q >= 2*k & k > 0)		// q >=k �� (log_3 k + 1)��O2,  
		{
			q = q - 2*k;				// 2->0
			k = k / 3;				// �ˬd�U�@�� 
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
		while (q >= 4 * k & k > 0)		// q >=k �� (log_5 k + 1)��O4,  
		{
			q = q - 4 * k;				// 4->0
			k = k / 5;				// �ˬd�U�@�� 
		}
		q = q + k;
	}
	return 0;
}

void FFT::swap(int &a,int &b)
{
	a ^= b;
	b ^= a;
	a ^= b;
}

int main()
{
	FFT t;
	//t.test();
	t.BitReverse5(125);
	cin.get();
	//system("pause");
	return 0;
}

