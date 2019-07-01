#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

void Gen_mat(const int n, double *a)
{
	srand(20190701);

//	#pragma omp parallel for
	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			a[i + j*n] = (double)rand() / RAND_MAX;
}

void Copy_mat(const int n, double *a, double *b)
{
//	#pragma omp parallel for
	#pragma omp for
	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			b[i + j*n] = a[i + j*n];
}

void Set_Iden(const int n, double *a)
{
//	#pragma omp parallel for
	#pragma omp for    
	for (int i=0; i<n; i++)
		for (int j=0; j<=i; j++)
			(i != j) ? a[i + j*n] = a[j + i*n] = 0.0: a[i + j*n] = 1.0;
}
