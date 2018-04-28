// Sorting_v2.cpp : 定義主控台應用程式的進入點。
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
	struct Node *prev;
};

class Sort{
public:
	void getdataL();
	void getdataA();
	void assign(Node*& cur, Node*& last, int a);
	void print_list();
	void print_array();
	void SelectSort();
private:
	int ncase, tempdata;
	int length=0;
	int *dataarray = NULL;
	Node *first, *end, *tmp;
	int iA = 0; //for array
	const int sizeMax = 1001;
};

void Sort::assign(Node*& cur, Node*& last, int a){
	cur->data = a;
	cur->prev = last;
	cur->next = NULL;
}

void Sort::getdataA(){
	ifstream ifs("Input.txt", ios::in);
	if (!ifs)
		cout << "Fail to open the input file." << endl;
	else
	{
		while (ifs >> tempdata)
			iA++;
	}
	ifs.close();
	ifstream ifs2("Input.txt", ios::in);
	if (!ifs2)
		cout << "Fail to open the input file." << endl;
	else
	{
		int c = 0;
		dataarray = new int[iA];
		while (!ifs2.eof()){
			ifs2 >> tempdata;
			if (tempdata != NULL){
				dataarray[c++] = tempdata;
			}
		}
	}
	ifs2.close();
}

void Sort::getdataL(){
	ifstream ifs("Input.txt", ios::in);
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
			if (first == NULL){
				Node *now = new Node;
				assign(now, tmp, tempdata);
				first = now;
				tmp = now;
			}
			else{
				Node *now = new Node;
				assign(now, tmp, tempdata);
				tmp->next = now;
				tmp = tmp->next;
			}
			ifs >> tempdata;
			//cout << tempdata <<" ";
		}
		Node *now = new Node;
		assign(now, tmp, tempdata);
		tmp->next = now;
		tmp = tmp->next;
		end = now;
	}
	ifs.close();
	if ((length & 1) == 0) //even
		ncase = 2;
	else
		ncase = 1;
}

void Sort::print_list(){
	tmp = first;
	do{
		ofs << tmp->data << "\t";
		cout << tmp->data << "\t";
		tmp = tmp->next;
	} while (tmp != NULL);
	cout << endl;
}

void Sort::print_array(){
	int c;
	for (c = 0; c < iA; c++){
		cout << dataarray[c] << "\t";
		ofs << dataarray[c] << "\t";
	}
	cout << endl;
	ofs << endl;
}

int main()
{
	Sort test;
	cout << "Array" << endl;
	test.getdataA();
	test.print_array();
	cout << "Linked list" << endl;
	test.getdataL();
	test.print_list();


	system("pause");
	return 0;
}

