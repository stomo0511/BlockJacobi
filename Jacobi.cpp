#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <mkl.h>
#include "Jacobi.hpp"

using namespace std;

///
/// @fn
/// Generate test matrix
/// @param (mode) mode (1 to 5)
/// @param (cond) condition number
/// @param (dmax) maximum avs. value of singular values
/// @param (n)    size of square matrix
/// @param (D)    pointer to the matrix
///
void Gen_test_mat(int mode, double cond, double dmax, const int n, double* A)
{
	assert( mode > 0 && mode < 6 );
	assert( cond > 0.0 );

	int iseed[4];
	double* sv = new double[n];

	iseed[0] = 2019; iseed[1] = 07; iseed[2] = 07; iseed[3] = 4093;

	LAPACKE_dlatms ( LAPACK_COL_MAJOR, n, n, 'U', iseed, 'N',
			sv, mode, cond, dmax, n-1, n-1, 'N', A, n);

	delete [] sv;
}

// Generate random matrix
void Gen_rand_mat(const int m, const int n, double *A)
{
	srand(20190611);

//	#pragma omp parallel for
	for (int i=0; i<m*n; i++)
		A[i] = 1.0 - 2*(double)rand() / RAND_MAX;
}

// Copy matrix A to B
void Copy_mat(const int m, const int n, double *A, double *B)
{
//	#pragma omp parallel for
	#pragma omp for
	for (int i=0; i<m*n; i++)
		B[i] = A[i];
}

// Show matrix
void Show_mat(const int m, const int n, double *A)
{
	cout.setf(ios::scientific);
	for (int i=0; i<m; i++) {
		for (int j=0; j<n; j++)
			cout << showpos << setprecision(4) << A[i + j*m] << ", ";
		cout << endl;
	}
	cout << endl;
}

// Set the matrix to identity matrix
void Set_Iden(const int n, double *A) // for square matrix only
{
//	#pragma omp parallel for
	#pragma omp for    
	for (int i=0; i<n; i++)
		for (int j=0; j<=i; j++)
			(i != j) ? A[i + j*n] = A[j + i*n] = 0.0: A[i + j*n] = 1.0;
}

// Norm of the off-diagonal elements
double Off_d(const int n, double *A) // for square matrix only
{
	double tmp = 0.0;

	#pragma omp parallel for reduction(+:tmp)
	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			if (j != i)
				tmp += A[i + j*n] * A[i + j*n];
	return sqrt(tmp);
}

void music(const int n, int *top, int *bot)
{
	int m = n/2;
	int *ct = new int[m];
	int *cb = new int[m];

	for (int i=0; i<m; i++)
	{
		ct[i] = top[i]; cb[i] = bot[i];
	}

	for (int k=0; k<m; k++)
	{
		if (k==0)
			top[0] = 0;
		else if (k==1)
			top[1] = cb[0];
		else if (k>1)
			top[k] = ct[k-1];
		if (k==m-1)
			bot[k] = ct[k];
		else
			bot[k] = cb[k+1];
	}
}

bool modulus_pair(const int w, const int i, const int r, int &p, int &q)
{
	assert(r <= w);
	assert(i <= w);

	if (r+1 < w - 2*i) {   // r <- r+1
		p = i+r; q = w-i-1;
		return false;
	} else if (r+1 == w - 2*i) {   // r <- r+1
		p = i+r - w/2; q = i+r;
		return true;
	} else {
		p = w/2 - i-1; q = i+r - w/2;   // w/2 -i <- w/2 -i-1
		return false;
	}
}
