// hw1_e14023056.cpp : 
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <windows.h>
using namespace std;
string ofname;
class Matrix
{
public:
	int **buildMatrix(int row, int column);
	int **getMatrix() { return matrix; }
	int multiplication();
	void printA();
private:
	int i, j;
	int size_Ar, size_Ac;
	int size_xr, size_xc;
	int **matrix;
};

int **Matrix::buildMatrix(int row, int column) {
	size_Ar = row;
	size_Ac = column;
	matrix = new int*[row];
	for (i = 0; i < row; i++)
	{
		matrix[i] = new int[column];
	}
	srand(time(NULL));
	#pragma omp parallel for
	for (i = 0; i < row; i++)
	{
	#pragma omp parallel for
		for (j = 0; j < column; j++)
		{
			matrix[i][j] = rand();
		}	
	}
	return matrix;
}
void Matrix::printA()
{
	ofstream ofs(ofname, ios::out);
	if (!ofs)
		ofs << "Fail to open the ouput file." << endl;
	else
	{
		for (i = 0; i < size_Ar; i++)
		{
			for (j = 0; j < size_Ac; j++)
			{
				cout << matrix[i][j] << "\t";
				ofs << matrix[i][j] << "\t";
			}
			cout << endl;
			ofs << endl;
		}
	}
}

int main()
{
	//clock_t tic, toc; //time_start and time_end
	double T0, T1,T2;
	int i, j;
	int row = 10, column = 10;
	int **matrix1, **matrix2;
	Matrix test;
	ofname = "Result.txt";
	matrix2 = test.buildMatrix(10,10);
	test.printA();
	
	/*
	tic = clock();
	for (i = 0; i < row; i++)
	{
	}
	toc = clock();
	T0 = (toc - tic) / (double)(CLOCKS_PER_SEC);
	
	tic = clock();
	matrix1 = new int*[row];
	for (i = 0; i < row; i++)
	{
		matrix1[i] = new int[column];
	}
	toc = clock();
	T1 = (toc - tic) / (double)(CLOCKS_PER_SEC);
	printf("settling time(wo/ omp) = \t%1f (s) \n", (T1 - T0) / 2.0);

	//settling data[row][column]
	tic = clock();
	matrix2 = new int*[row];
	#pragma omp parallel for
	for (i = 0; i < row; i++)
	{
		matrix2[i] = new int[column];
	}
	toc = clock();
	T1 = (toc - tic) / (double)(CLOCKS_PER_SEC);
	printf("settling time(w/ omp) = \t%1f (s) \n\n", (T1 - T0) / 2.0);

	srand(time(NULL));

	
	tic = clock();
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < column; j++)
		{
			matrix1[i][j] = rand();
			cout << matrix1[i][j] << "\t";
		}
		cout << endl;
	}
	toc = clock();
	T2 = (toc - tic) / (double)(CLOCKS_PER_SEC);
	printf("assignment time(wo/ omp) = \t%1f (s) \n", (T2 - T0) / 2.0);
	
	//data assignment
	tic = clock();
	#pragma omp parallel for
		for (i = 0; i < row; i++)
		{
			#pragma omp parallel for
			for (j = 0; j < column; j++)
			{
				matrix2[i][j] = rand();
				cout << matrix2[i][j] << "\t";
			}
			cout << endl;
		}
		toc = clock();
		T2 = (toc - tic) / (double)(CLOCKS_PER_SEC);
		printf("assignment time(w/ omp) = \t%1f (s) \n\n", (T2 - T0) / 2.0);

	//release
		for (i = 0; i < row; i++)
		{
			delete[] matrix1[i];
			delete[] matrix2[i];
		}
	delete[] matrix1;
	delete[] matrix2;
	*/
	system("pause");
	return 0;
}