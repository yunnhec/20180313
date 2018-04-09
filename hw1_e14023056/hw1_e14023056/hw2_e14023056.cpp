// hw2_e14023056.cpp : 
//

#include "stdafx.h"
#include <iostream>
#include <time.h>
#include <math.h>
#include <windows.h>
using namespace std;

int main()
{
	clock_t tic, toc; //time_start and time_end
	int i, j;
	int row = 10, column = 10;
	int **matrix;

	tic = clock();
	matrix = new int*[row];
	for (i = 0; i < row; i++)
	{
		matrix[i] = new int[column];
	}
	toc = clock();
	T0 = (toc - tic) / (double)(CLOCKS_PER_SEC);

	//settling data[row][column]
	tic = clock();
	matrix = new int*[row];
	for (i = 0; i < row; i++)
	{
		matrix[i] = new int[column];
	}
	toc = clock();
	T1 = (toc - tic) / (double)(CLOCKS_PER_SEC);
	printf("settling time = \t%1f (s) \n\n", (T1 - T0) / 2.0);

	
	srand(time(NULL));

	//data assignment
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < column; j++)
		{
			matrix[i][j] = rand();
			//cout << matrix[i][j] << "\t";
		}
		//cout << endl;
	}

	//release
	for (i = 0; i < row; i++)
		delete[] matrix[i];
	delete[] matrix;

	system("pause");
	return 0;
}