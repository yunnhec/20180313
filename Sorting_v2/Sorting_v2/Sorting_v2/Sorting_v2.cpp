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
	void getdata();
	void assign(Node*& cur, Node*& last, int a);
	void print_list();
private:
	int ncase, tempdata;
	int length=0;
	Node *first, *end, *tmp;
};

void Sort::assign(Node*& cur, Node*& last, int a){
	cur->data = a;
	cur->prev = last;
	cur->next = NULL;
}

void Sort::getdata(){
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

void Sort::print_list(){
	tmp = first;
	do{
		//ofs << tmp->data << " ";
		cout << tmp->data << " ";
		tmp = tmp->next;
	} while (tmp != NULL);
	cout << endl;
}

int main()
{
	Sort test;
	test.getdata();
	test.print_list();

	system("pause");
	return 0;
}

