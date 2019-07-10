//
// Blocked version of One sided Jacobi
//

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <mkl_cblas.h>
#include "Jacobi.hpp"

using namespace std;

#define CASHSIZE 24000

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
	Gen_test_mat( 3, 100001.0, n, G);
//	Show_mat(n,n,G);

//	cout << noshowpos << "p = " << p << showpos << endl;
//	for (int i=0; i<n; i++) {
//		for (int j=0; j<nb; j++)
//			cout << G[i + (p*nb+j)*n] << ", ";
//		cout << endl;
//	}
//	cout << endl << noshowpos << "q = " << q << showpos << endl;
//	for (int i=0; i<n; i++) {
//		for (int j=0; j<nb; j++)
//			cout << G[i + (q*nb+j)*n] << ", ";
//		cout << endl;
//	}
//	cout << endl;

	const int nb1 = CASHSIZE / n;
	assert(CASHSIZE % n == 0);

//	#pragma omp parallel
//	{
//		const int nt = omp_get_num_threads();
//		const int id = omp_get_thread_num();
//		const int nb = n/nt;
//		assert(n % nt == 0);
//
//		int pp, qq;
//		bool flg;
//
//		for (int r=0; r<2*nt; r++)
//		{
//			flg = modulus_pair(2*nt,id,r,pp,qq);
//
////			#pragma omp critical
////			cout << "id = " << id << ": r = " << r << ", (pp,qq) = (" << pp << ", " << qq << ")\n";
//
//			if (flg == false)
//			{
//				for (int p2=nb*(pp-1); p2<nb*pp; p2+=nb1)
//					for (int q=nb*(qq-1); q<nb*qq; q++)
//						for (int p=p2; p<p2+nb1-1; p++)
//							cout << "(p,q) = (" << p << ", " << q << ")\n";
//			}
//			else
//			{
//				for (int p2=nb*(pp-1)+1; p2<nb*pp; p2+=nb1)
//					for (int q=p2+1; q<nb*pp; q++)
//						for (int p=p2; p<q; p++)
//							cout << "(p,q) = (" << p << ", " << q << ")\n";
//				for (int p2=nb*(pp-1)+1; p2<nb*qq; p2+=nb1)
//					for (int q=p2+1; q<nb*pp; q++)
//						for (int p=p2; p<q; p++)
//							cout << "(p,q) = (" << p << ", " << q << ")\n";
//			}
//		}
//	}


	delete [] G;

	return EXIT_SUCCESS;
}
