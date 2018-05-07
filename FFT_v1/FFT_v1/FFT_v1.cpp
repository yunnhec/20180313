// FFT_v1.cpp : 定義主控台應用程式的進入點。
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
	int BitReverse(int num);
	void swap(int &a, int &b);
	void test();
};

void FFT::test(){
	int x = 128, y = 1;
	cout << "(x,y) = " << x << ", " << y << endl;
	swap(x, y);
	cout << "(x,y) = " << x << ", " << y << endl;
}

int FFT::BitReverse(int num){

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
	t.test();
	cin.get();
	//system("pause");
	return 0;
}

