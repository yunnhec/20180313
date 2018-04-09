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
	void initialize();	//pseudo-random number generator
	int **buildMatrixA(int rowA, int columnA); //genetare a matrix with size (rowA, columnA)
	int **buildMatrixX(int columnX); //genetare a matrix with size (columnA, columnX)
	int **getMatrixA() { return matrixA; } //return matrixA
	int **getMatrixX() { return matrixX; } //return matrixX
	int **Multiplication(); //implement matrix multiplication and return a matrix of A*X
	void print(); //print all result to screen and file
	void printfile(); //print all result to file
	void release(); //release memory
private:
	int i, j, k;
	int size_Ar, size_Ac;
	int size_Xr, size_Xc;
	int **matrixA, **matrixX, **matrixB;
};

void Matrix::initialize()
{
	srand(time(NULL));
}

int **Matrix::buildMatrixA(int rowA, int columnA) {
	size_Ar = rowA;
	size_Ac = columnA;

	matrixA = new int*[size_Ar];
	#pragma omp parallel for
	for (i = 0; i < size_Ar; i++)
	{
		matrixA[i] = new int[size_Ac];
	}
	
	#pragma omp parallel for
	for (i = 0; i < size_Ar; i++)
	{
	#pragma omp parallel for
		for (j = 0; j < size_Ac; j++)
		{
			matrixA[i][j] = rand()%10; //all values are in 0~10
		}	
	}
	return matrixA;
}

int **Matrix::buildMatrixX(int columnX) {
	int rowX = size_Ac;
	size_Xr = rowX;
	size_Xc = columnX;

	matrixX = new int*[size_Xr];
	#pragma omp parallel for
	for (i = 0; i < size_Xr; i++)
	{
		matrixX[i] = new int[size_Xc];
	}
	#pragma omp parallel for
	for (i = 0; i < size_Xr; i++)
	{
		#pragma omp parallel for
		for (j = 0; j < size_Xc; j++)
		{
			matrixX[i][j] = rand()%10; //all values are in 0~10
		}
	}
	return matrixX;
}

int **Matrix::Multiplication()
{
	matrixB = new int*[size_Ar];
	#pragma omp parallel for
	for (i = 0; i < size_Ar; i++)
	{
		matrixB[i] = new int[size_Xc];
	}
	#pragma omp parallel shared(matrixA,matrixB,matrixX) private(i,j,k)  
	{
		#pragma omp for schedule(dynamic)  
		for (i = 0; i < size_Ar; i++)
		{
			for (j = 0; j < size_Xc; j++)
			{
				matrixB[i][j] = 0;
				for (k = 0; k < size_Ac; k++)
				{
					matrixB[i][j] += matrixA[i][k] * matrixX[k][j];
				}
			}
		}
	}
	return matrixB;
}

void Matrix::print()
{
	ofstream ofs(ofname, ios::out);
	if (!ofs)
		ofs << "Fail to open the ouput file." << endl;
	else
	{
		cout << "matrix A =" << endl;
		ofs << "A=[";
		for (i = 0; i < size_Ar; i++)
		{
			for (j = 0; j < size_Ac; j++)
			{
				cout << matrixA[i][j] << "\t";
				ofs << matrixA[i][j] << "\t";
			}
			cout << endl;
			ofs<< endl;
		}
		cout  << endl << "matrix X =" << endl;
		ofs << "];" << endl << "X=[";
		for (i = 0; i < size_Xr; i++)
		{
			for (j = 0; j < size_Xc; j++)
			{
				cout << matrixX[i][j] << "\t";
				ofs << matrixX[i][j] << "\t";
			}
			cout << endl;
			ofs << endl;
		}
		cout << endl << "A * X =" << endl;
		ofs << "];" << endl << "A * X =" << endl;
		for (i = 0; i < size_Ar; i++)
		{
			for (j = 0; j < size_Xc; j++)
			{
				cout << matrixB[i][j] << "\t";
				ofs << matrixB[i][j] << "\t";
			}
			cout << endl;
			ofs << endl;
		}
	}
}

void Matrix::printfile()
{
	ofstream ofs(ofname, ios::out);
	if (!ofs)
		ofs << "Fail to open the ouput file." << endl;
	else
	{
		ofs << "A=[";
		for (i = 0; i < size_Ar; i++)
		{
			for (j = 0; j < size_Ac; j++)
			{
				ofs << matrixA[i][j] << "\t";
			}
			ofs << endl;
		}
		ofs << "];" << endl << "X=[";
		for (i = 0; i < size_Xr; i++)
		{
			for (j = 0; j < size_Xc; j++)
			{
				ofs << matrixX[i][j] << "\t";
			}
			ofs << endl;
		}
		ofs << "];" << endl << "A * X =" << endl;
		for (i = 0; i < size_Ar; i++)
		{
			for (j = 0; j < size_Xc; j++)
			{
				ofs << matrixB[i][j] << "\t";
			}
			ofs << endl;
		}
	}
}

void Matrix::release()
{
	for (i = 0; i < size_Ar; i++)
	{
		delete[] matrixA[i];
	}
	delete[] matrixA;
	for (i = 0; i < size_Xr; i++)
	{
		delete[] matrixX[i];
	}
	delete[] matrixX;
	for (i = 0; i < size_Ar; i++)
	{
		delete[] matrixB[i];
	}
	delete[] matrixB;
}

int main()
{
	clock_t tic, toc; //time_start and time_end
	double T0, T1;
	int i, j,k;
	int row = 1000, column = 1000, column2 = 1;
	int **A, **X, **B;
	Matrix test;
	ofname = "Result.txt";

	tic = clock();
	#pragma omp parallel private(i,j,k)  
	{
	#pragma omp for schedule(dynamic)  
		for (i = 0; i < row; i++){
			for (j = 0; j < column2; j++){
				for (k = 0; k < column; k++){
				}
			}
		}
	}
	toc = clock();
	T0 = (toc - tic) / (double)(CLOCKS_PER_SEC);

	test.initialize();
	A = test.buildMatrixA(row,column);
	X = test.buildMatrixX(column2);

	tic = clock();
	B = test.Multiplication();
	toc = clock();

	test.printfile();

	T1 = (toc - tic) / (double)(CLOCKS_PER_SEC);
	printf("Operation time = %f\n" ,T1-T0);
	test.release();

	system("pause");
	return 0;
}