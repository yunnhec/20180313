// Sorting.cpp : Starting of the program
//

#include "stdafx.h"
#include <Windows.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;
ofstream ofs("Output.txt", ios::out);

struct Node
{
	int data;
	struct Node *next;
};

class QuickSort
{
public:
	void ReadData();
	void assign(Node*& cur, int a);
	void swap(Node*& a, Node*& b);
	void print();
	void QuickMedian();
private:
	int pivot;
	int tempdata;
	Node* tmp;
	Node* list; //first
	Node* now;
	int length=0; //number of data
};

void QuickSort::assign(Node*& cur, int a)
{
	cur->data = a;
	cur->next = NULL;
}

void QuickSort::ReadData()
{
	//list = NULL; //initialize
	ifstream ifs("Input1.txt", ios::in);
	if (!ifs)
		cout << "Fail to open the input file." << endl;
	else
	{
		//first element
		ifs >> tempdata;
		//cout << tempdata << " ";
		length++;
		while (!ifs.eof()){
			if (tempdata == '\n')
				break;
			length++;
			if (list == NULL){
				Node *now = new Node;
				assign(now, tempdata);
				list = now;
				tmp = now;
			}
			else{
				Node *now = new Node;
				assign(now, tempdata);
				tmp->next = now;
				tmp = tmp->next;
			}
			ifs >> tempdata;
			//cout << tempdata <<" ";
		}
		Node *now = new Node;
		assign(now, tempdata);
		tmp->next = now;
		tmp = tmp->next;
	}
}

void QuickSort::swap(Node*& m, Node*& n)
{
	Node* swaptmp = m;
	m = n;
	n = swaptmp;
}

void QuickSort::print()
{
	now = list;
	do{
		ofs << now->data << " ";
		cout << now->data << " ";
		now = now->next;
	} while (now != NULL);
	cout << endl;
}

int main()
{
	QuickSort test;
	test.ReadData();
	test.print();
	system("pause");
	return 0;
}


