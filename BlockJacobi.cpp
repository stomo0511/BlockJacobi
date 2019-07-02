#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <mkl_cblas.h>
#include "Jacobi.hpp"

using namespace std;

int main(int argc, char **argv)
{
	// Usage "a.out [size of matrix] [size of submatrix]"
	assert(argc > 2);

	// Matrix size: n x n
	const int n = atoi(argv[1]);
	assert(n % 2 == 0);

	// Sub-matrix size: nb x nb
	const int nb = atoi(argv[2]);
	assert(n % nb == 0);

	// # of sub-matrix
	const int w = n / nb;
	assert(w % 2 == 0);

	double *G = new double[n*n];   // Original matrix

	Gen_mat(n,n,G);
	Show_mat(n,n,G);

	int p = 0;
	int q = 2;

	cout << noshowpos << "p = " << p << showpos << endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<nb; j++)
			cout << G[i + (p*nb+j)*n] << ", ";
		cout << endl;
	}
	cout << endl << noshowpos << "q = " << q << showpos << endl;
	for (int i=0; i<n; i++) {
		for (int j=0; j<nb; j++)
			cout << G[i + (q*nb+j)*n] << ", ";
		cout << endl;
	}
	cout << endl;



	delete [] G;

	return EXIT_SUCCESS;
}
