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
	int MedianA();
	void swapA(int s1, int s2);
	int sortA(int *p, int range_L, int range_R);
private:
	int tempdata, result = 0;
	int length=0;
	int *dataarray = NULL;
	Node *first, *end, *tmp, *med;
	int iA = 0; //for array
	const int sizeMax = 1001;
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
	ifstream ifs("Input1.txt", ios::in);
	if (!ifs)
		cout << "Fail to open the input file." << endl;
	else
	{
		while (ifs >> tempdata)
			iA++;
	}
	ifs.close();
	ifstream ifs2("Input1.txt", ios::in);
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
	int pivot = SelectionA(range_L, (range_R - range_L) / 2, range_R);
	cout << "pivot_index = " << pivot << endl;
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
		cout << "left = " << left << " ,right = " << right << endl;
		swapA(left, right);
		print_array();
	}
	if (selCase){
		swapA(pivot, left);
		pivot = left;
	}
	else{
		swapA(pivot, right);
		pivot = right;
	}
	swapA(pivot, left);
	cout << "left = " << left << " ,right = " << right << endl;
	print_array();
	cout << endl;
	//先只考慮奇數
	if (pivot == iA / 2)
		return pivot;
	else if (pivot > iA / 2) {  //pivot在數列右邊->取左半
		sortA(dataarray, tmpL, pivot - 1);
	}
	else //pivot在數列左邊->取右半
	{
		sortA(dataarray, pivot + 1, tmpR);
	}
	//return pivot;
	return 0;
}

int Sort::MedianA(){
	int p = sortA(dataarray, 0, iA - 1);
	//確認pivot的位置
	if (p == iA / 2){
		if ((iA & 1) == 1){ //odd
			result = p;
			return result;
		}
		else{ //even
			int p2 = sortA(dataarray, 0, iA - 1);
		}
	}
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
	cout << "Array" << endl;
	test.getdataA();
	test.print_array();
	cout << "Linked list" << endl;
	test.getdataL();
	test.print_list();
	cout << endl;
	test.MedianA();
	//cout << "Median = " << test.MedianA() << endl;
	system("pause");
	return 0;
}

