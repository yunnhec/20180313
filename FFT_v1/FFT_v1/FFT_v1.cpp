// FFT_v1.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include "fft_v1.h"

#include <time.h>
#include <omp.h>
using namespace std;

int main()
{
	clock_t t1, t2;
	FFT t;
	t1 = clock();
	t.fft(2, 0, 0);
	t2 = clock();
	
	printf("time = %f\n", (t2 - t1) / (double)(CLOCKS_PER_SEC));
	system("pause");
	t.getX();
	system("pause");
	return 0;
}

