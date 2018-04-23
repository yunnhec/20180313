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
	struct Node *prev;
};

class QuickSort
{
public:
	void ReadData();
	void assign(Node*& cur, Node*& last, int a);
	void swap(Node*& a, Node*& b);
	void print();
	void sort(Node*& p, int nnr, int nnl);
	void QuickMedian();
	int result();
private:
	int nr, nl;
	int med=0;
	Node* pivot, *left, *right;
	int tempdata;
	Node* tmp;
	Node* first, *end;
	Node* now;
	int length=0,ncase=0; //number of data
};

void QuickSort::assign(Node*& cur,Node*& last, int a)
{
	cur->data = a;
	cur->prev = last;
	cur->next = NULL;
}

void QuickSort::ReadData()
{
	//first = NULL; //initialize
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
	if ((length & 1) == 0) //even
		ncase = 2;
	else
		ncase = 1;
}

void QuickSort::sort(Node*& p, int nnr, int nnl){
	if (first != pivot)
		left = first;
	else
	{
		left = first->next;
		nnl++;
	}
	right = end;
	int tmpnr = nnr;
	int tmpnl = nnl;
	//find left > pivot
	while (!(left->data > p->data)){ 
		if (left->next != NULL){
			left = left->next;
			nnl++;
		}
		else
			nnl = tmpnl;	
	}
	//find right <= pivot
	while (!(right->data <= p->data)){
		if (right->prev != NULL){
		right = right->prev;
		nnr--;
		}
		else
			nnr = tmpnr;
	}
	if (nnr >= nnl)
		swap(right, left);
	else{
		swap(right, p);
		nr = nnr;
		nl = nnl;
	}
	cout << "nnr = " << nnr;
	cout << ", nnl = " << nnl << endl;
	cout << endl;
}

void QuickSort::QuickMedian(){
	pivot = first;
	nr = length;
	nl = 1;
	cout << "nr = " << nr;
	cout << ", nl = " << nl <<", p = ";
	cout << pivot->data << endl;
	cout << endl;
	cout << endl;
	sort(pivot, nr, nl);
	cout << "nr = " << nr;
	cout << ", nl = " << nl << ", p = ";
	cout << pivot->data << endl;
	cout << endl;
	print();
	//-------------------------------------------
	while (nl == length || nr ==1)
	{
		if (nr > length/2){
			
			sort(pivot, nr, nl); //nl
		}
		else{
			sort(pivot, nr, nl);
		}
		//-------------------------------------------
		//cout << "nr = " << nr << endl;
		print();
		switch (ncase){
		case '1':
			if (nr == (length / 2 + 1)) {//odd
				med = pivot->data;
				break;
			}
		case '2':
			if (nr == length / 2) {//even
				tmp = pivot->prev;
				med = (pivot->data + tmp->data) / 2;
				break;
			}
		default:
			continue;
		}
		break;
	}
}

void QuickSort::swap(Node*& m, Node*& n)
{
	int swaptmp = m->data;
	m->data = n->data;
	n->data = swaptmp;
}

void QuickSort::print()
{
	now = first;
	do{
		ofs << now->data << " ";
		cout << now->data << " ";
		now = now->next;
	} while (now != NULL);
	cout << endl;
	/*
	now = end;
	do{
		ofs << now->data << " ";
		cout << now->data << " ";
		now = now->prev;
	} while (now != NULL);
	cout << endl;
	ofs << endl;*/
}

int QuickSort::result(){
	/*if (nr == (length / 2 + 1)) //odd
		med = pivot->data;
	else if (nr == length / 2) {//even
		tmp = pivot->prev;
		med = (pivot->data + tmp->data) / 2;
	}*/
	return med;
}

int main()
{
	QuickSort test;
	test.ReadData();
	test.print();
	cout << endl;
	cout << endl;
	cout << endl;
	test.QuickMedian();
	test.print();
	cout << test.result() << endl;
	system("pause");
	return 0;
}


