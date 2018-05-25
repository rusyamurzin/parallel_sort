#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "fstream"
#include "omp.h"
#include "pthread.h"
#include <ctime>
#include "Math.h"
using namespace std;
int ThreadNum; // Number of threads
int swapped = 0;
int *c;
int dim;

void BubbleSort(int *A, int start, int end) {
	for (int i = start; i < end; i++){
		for (int j = start; j < start + end - i; j++) {
			if (A[j] > A[j+1]) {
				int temp = A[j];
				A[j] = A[j+1];
				A[j+1] = temp;
			}
		}
	}
}
// Function for merging of two sorted blocks
void MergeBlocks(int* pData, int Index1, int BlockSize1, int Index2, int BlockSize2) {
	int* pTempArray = new int[BlockSize1 + BlockSize2];
	int i1 = Index1, i2 = Index2, curr = 0;
	while ((i1 < (Index1 + BlockSize1)) && (i2 < (Index2 + BlockSize2))) {
		if (pData[i1] < pData[i2]) 
			pTempArray[curr++] = pData[i1++];
		else {
			pTempArray[curr++] = pData[i2++];
			swapped = 1;
		}
	}
	while (i1 < (Index1 + BlockSize1))
		pTempArray[curr++] = pData[i1++];
	while (i2 < (Index2 + BlockSize2))
		pTempArray[curr++] = pData[i2++];
	for (int i = 0; i < BlockSize1 + BlockSize2; i++)
		pData[Index1 + i] = pTempArray[i];
	delete[] pTempArray;
}
void* sortArr(void *num)
{
	int first = 0, last;
	int cast = (int)num;
	//Calculating first and last element
	if (ThreadNum == 0) {
		first = 0;
		last = dim - 1;
	}
	else {
		first = (cast * (dim / ThreadNum));
		last = first + (dim / ThreadNum) - 1;
	}

	BubbleSort(c, first, last);
	return NULL;
}
// Function for parallel odd-even transposotion
void PthreadOddEvenSort(int* pData, int Size, int num_thr) {
	pthread_t *thread = new pthread_t[num_thr];
	for (int i = 0; i < num_thr; i++) {
		pthread_create(&thread[i], NULL, sortArr, (void*)i);
	}
	for (int i = 0; i < num_thr; i++) {
		pthread_join(thread[i], NULL);
	}

	int* Index = new int[num_thr];
	int* BlockSize = new int[num_thr];
	for (int i = 0; i < num_thr; i++) {
		Index[i] = int((i*Size) / double(num_thr));
		if (i < num_thr - 1)
			BlockSize[i] = int(((i + 1)*Size) / double(num_thr)) - Index[i];
		else
			BlockSize[i] = Size - Index[i];
	}
	// Odd-even transposition of data blocks
	int Iter = 0;
	do {
		swapped = 0;
		int powIter = (int)pow(2, Iter);
		if(powIter < num_thr){
			for (int i = 0; i < num_thr; i += 2 * powIter) {
				int sizeOfLeftBlocks = 0;
				int sizeOfRightBlocks = 0;
				for (int j = 0; j < powIter; j++) {
					sizeOfLeftBlocks += BlockSize[i + j];
					sizeOfRightBlocks += BlockSize[i + powIter + j];
				}
				MergeBlocks(pData, Index[i], sizeOfLeftBlocks,
					Index[i + powIter], sizeOfRightBlocks);
			}
		}
		Iter++;
	} while (swapped);
	delete[] thread;
	delete[] Index;
	delete[] BlockSize;
}


// Function for checking if the array is sorted
bool IsSorted(int* pData, int Size) {
	bool res = true;
	for (int i = 1; (i<Size) && (res); i++) {
		if (pData[i]<pData[i - 1])
			res = false;
	}
	return res;
}


// Function for serial odd-even transposition
void OddEvenSort(int *pData, int Size) {
	int temp;
	int upper_bound;
	if (Size % 2 == 0)
		upper_bound = Size / 2 - 1;
	else
		upper_bound = Size / 2;
	for (int i = 0; i<Size; i++) {
		if (i % 2 == 0) // even iteration
			for (int j = 0; j < Size / 2; j++) {
				if (pData[2 * j] > pData[2 * j + 1]) {
					temp = pData[2 * j];
					pData[2 * j] = pData[2 * j + 1];
					pData[2 * j + 1] = temp;
				}
			}
		else // odd iteration
			for (int j = 0; j < upper_bound; j++) {
				if (pData[2 * j + 1] > pData[2 * j + 2]) {
					temp = pData[2 * j + 1];
					pData[2 * j + 1] = pData[2 * j + 2];
					pData[2 * j + 2] = temp;
				}
			}
	}
}
// Function for parallel odd-even transposotion
void OpenMPOddEvenSort(int* pData, int Size) {
	int upper_bound;
	if (Size % 2 == 0)
		upper_bound = Size / 2 - 1;
	else
		upper_bound = Size / 2;
	for (int i = 0; i<Size; i++) {
		if (i % 2 == 0) // Even iteration
#pragma omp parallel for
			for (int j = 0; j<Size / 2; j++){
				if (pData[2 * j] > pData[2 * j + 1]) {
					int temp = pData[2 * j];
					pData[2 * j] = pData[2 * j + 1];
					pData[2 * j + 1] = temp;
				}
			}
		else // Odd iteration
#pragma omp parallel for
			for (int j = 0; j < upper_bound; j++) {
				if (pData[2 * j + 1] > pData[2 * j + 2]) {
					int temp = pData[2 * j + 1];
					pData[2 * j + 1] = pData[2 * j + 2];
					pData[2 * j + 2] = temp;
				}
			}
	}
}


int main(int argc, char **argv)
{
	if (argc < 2) return 2;
	ifstream in(argv[1]);
	//ifstream in("D:\\parallel_sort\\x64\\Debug\\input50k.txt");
	if (!in.is_open()) cout << "Cannot open file" << endl;
	else cout << "File opened" << endl;
	int count = 0;
	in >> count;
	dim = count;
	cout << "Size: " << count << endl;
	int num_threads = 1;
	in >> num_threads;
	num_threads = 16;
	omp_set_num_threads(num_threads);
	ThreadNum = num_threads;
	cout << "Num threads: " << num_threads << endl;
	int *a = new int[count];
	for (int k = 0; k < count; k++)
		in >> a[k];
	int *b = new int[count];
	for (int k = 0; k < count; k++)
		b[k] = a[k];
	c = new int[count];
	for (int k = 0; k < count; k++)
		c[k] = a[k];

	//Serial OddEvenSort
	unsigned int start_time = clock();
	OddEvenSort(a, count);
	unsigned int end_time = clock();
	unsigned int search_time = end_time - start_time;
	printf("\n");
	cout << "Serial OddEvenSort time : " << search_time << " ms" << endl;

	//OpenMP OddEvenSort time
	start_time = clock();
	OpenMPOddEvenSort(b, count);
	end_time = clock();
	unsigned int search_OMP_time = end_time - start_time;
	cout << "OpenMP OddEvenSort time : " << search_OMP_time << " ms" << endl;

	//Pthread OddEvenSort
	start_time = clock();
	PthreadOddEvenSort(c, count, num_threads);
	end_time = clock();
	unsigned int search_Pthread_time = end_time - start_time;
	cout << "Pthread OddEvenSort time : " << search_Pthread_time << " ms" << endl;

	//for (int i = 0; i < 40; i++)
	//	printf("%d ", c[i]);
	printf("\n");
	/*
	std::ofstream vmdelet_out("output.txt", ios::app);
	vmdelet_out << search_time << endl;
	vmdelet_out.close();
	std::ofstream vmdelet_out1("output1.txt", ios::app);
	vmdelet_out1 << search_OMP_time << endl;
	vmdelet_out1.close();
	std::ofstream vmdelet_out2("output2.txt", ios::app);
	vmdelet_out2 << search_Pthread_time << endl;
	vmdelet_out2.close();
*/
	delete[] a;
	delete[] b;
	delete[] c;
	in.close();
	system("pause");

	return 0;
}
