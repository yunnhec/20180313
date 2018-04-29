// Sorting_v2.cpp : 定義主控台應用程式的進入點。
//

#include "stdafx.h"
#include <Windows.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
	Node* SelectionL(Node*& n1, Node*& n2, Node*& n3);
	Node *FindI(Node*& d, int index);
	int SelectionA(int n1, int n2, int n3);
	double MedianA();
	void swapA(int s1, int s2);
	int sortA(int *p, int range_L, int range_R);
private:
	int tempdata;
	double result = NULL;
	int length=0;
	int *dataarray = NULL;
	Node *first, *end, *tmp, *med;
	int iA = 0; //for array
};

void Sort::assign(Node*& cur, Node*& last, int a){
	cur->data = a;
	cur->prev = last;
	cur->next = NULL;
}

Node* Sort::FindI(Node*& d,int index){
	int c = 0;
	while (c++ != index)
		d = d->next;
	return d;
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
}

Node* Sort::SelectionL(Node*& n1, Node*& n2, Node*& n3){
	//return pivot Node
	int x1 = n1->data;
	int x2 = n2->data;
	int x3 = n3->data;
	if (x1 <= x2){
		if (x2 <= x3)
			return n3;
		else if (x3 <= x1)
			return n1;
		else
			return n2;
	}
	else{
		if (x1 <= x3)
			return n1;
		else if (x3 <= x2)
			return n2;
		else
			return n3;
	}
}

int  Sort::SelectionA(int n1, int n2, int n3){
	//return pivot index
	int x1 = dataarray[n1];
	int x2 = dataarray[n2];
	int x3 = dataarray[n3];
	if (x1 <= x2){
		if (x2 <= x3)
			return n2;
		else if (x3 <= x1)
			return n1;
		else
			return n3;
	}
	else{
		if (x1 <= x3)
			return n1;
		else if (x3 <= x2)
			return n2;
		else
			return n3;
	}
}

void Sort::swapA(int s1, int s2){
	int t = dataarray[s1];
	dataarray[s1] = dataarray[s2];
	dataarray[s2] = t;
}

int Sort::sortA(int *p, int range_L, int range_R){
	//choose a pivot
	int pivot = SelectionA(range_L, range_L+(range_R - range_L) / 2, range_R);
	int left, right,selCase;
	if (pivot == range_L){
		left = range_L+1;
		right = range_R;
		selCase = 0; //False
	}
	else if (pivot == range_R){
		left = range_L;
		right = range_R-1;
		selCase = 1; //True
	}
	else{
		left = range_L;
		right = range_R;
		selCase = 1; //True
	}
	int tmpL = left, tmpR = right;
	while (1) {
		left = tmpL;
		right = tmpR;
		//從左邊找>pivot的位置
		while (left < range_R && dataarray[pivot] > dataarray[left])
			left++;
		//從右邊找<=pivot的位置
		while (right>range_L && dataarray[pivot] <= dataarray[right])
			right--;
		if (left >= right) break;
		swapA(left, right);
	}
	if (selCase){
		swapA(pivot, left);
		pivot = left;
	}
	else{
		swapA(pivot, right);
		pivot = right;
	}
	if ((iA & 1) == 1){ //odd
		if (pivot == (iA) / 2)
			return pivot;
		else if (pivot > (iA) / 2) //pivot在數列右邊->取左半
			sortA(dataarray, tmpL, pivot - 1);
		else //pivot在數列左邊->取右半
			sortA(dataarray, pivot + 1, tmpR);
	}
	else{//even
		if (pivot == (iA-1) / 2)
			return pivot;
		else if (pivot > (iA-1) / 2) //pivot在數列右邊->取左半
			sortA(dataarray, tmpL, pivot - 1);
		else //pivot在數列左邊->取右半
			sortA(dataarray, pivot + 1, tmpR);
	}
}

double Sort::MedianA(){
	int p1 = sortA(dataarray, 0, iA - 1);
	if ((iA & 1) == 0){ //even
		int smallest=dataarray[p1+1];
		for (int c = p1 + 1; c < iA; c++){
			if (dataarray[c] < smallest)
				smallest = dataarray[c];
		}
		result = (dataarray[p1]+smallest)/2.0;
	}
	else{ //odd
		result = (double)dataarray[p1];
	}
	return result;
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

void Sort::print_list(){
	tmp = first;
	do{
		ofs << tmp->data << "\t";
		cout << tmp->data << "\t";
		tmp = tmp->next;
	} while (tmp != NULL);
	cout << endl;
}

int main()
{
	Sort test;
	cout << "Input data: " << endl;
	ofs << "Input data: " << endl;
	test.getdataA();
	test.print_array();
	cout << endl;
	ofs << endl;
	test.MedianA();
	cout << "Median = " << test.MedianA() << endl;
	ofs << "Median = " << test.MedianA() << endl;
	system("pause");
	return 0;
}

