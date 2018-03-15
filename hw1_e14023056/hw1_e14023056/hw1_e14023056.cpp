// hw1_e14023056.cpp : 
//

#include "stdafx.h"
#include <time.h>
#include <math.h>
#include <windows.h>
using namespace std;

int main()
{
	clock_t tic, toc; //time_start and time_end
	int N = 1000000000; // 1e9
	double T0, T1, T2, T3;

	tic = clock();
	double a = 3.547564, b = 7.55448, c = 1.235484;

	printf("Running for %d times.\n\n", N);

	// time for running an empty loop
	for (int i = 0; i < N; i++)
	{
	}
	toc = clock();
	T0 = (toc - tic) / (double)(CLOCKS_PER_SEC);

	// time for running addition
	tic = clock();
	for (int i = 0; i < N; i++)
	{
		a = b + a;
		b = a - b;
		a = a - b;
	}
	toc = clock();
	T1 = (toc - tic) / (double)(CLOCKS_PER_SEC);
	printf("1. T(addition) = \t%1f (s) \n", (T1 - T0)/2.0 );
	printf("(a,b) = (%f , %f) \n\n", a, b);

	// time for running multiplication
	a = 6.547564, b = 1.55448;

	for (int i = 0; i < N; i++)
	{
		a = b / c;
		b = a / b;
		c = a * b;
	}
	toc = clock();
	T2 = (toc - tic) / (double)(CLOCKS_PER_SEC);
	printf("2. T(multiplication) = \t %1f (s) \n", (T2 - T0)/3.0 );
	printf("(a,b) = (%f , %f) \n\n", a, b);

	// time for calculating sin
	a = 2.547564, b = 4.55448;
	tic = clock();
	for (int i = 0; i < N; i++)
	{
		a = sin(b);
		b = sin(a);
	}
	toc = clock();
	T3 = (toc - tic) / (double)(CLOCKS_PER_SEC);
	printf("3. T(sin) = \t %1f (s) \n", (T3 - T0)/2.0 );
	printf("(a,b) = (%f , %f) \n\n", a, b);

	system("pause");
	return 0;
}

