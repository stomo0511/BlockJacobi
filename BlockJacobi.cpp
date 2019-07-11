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

	double *A = new double[n*n];   // Original matrix
	double *V = new double[n*n];   // Right-singular vector matrix

	int *top = new int[n/2];
	int *bot = new int[n/2];

	Gen_test_mat( 5, 10000.0, 100.0, n, A);

	#pragma omp parallel
	{
		Set_Iden(n,V);     // V <- I

		#pragma omp for
		for (int i=0;i<n/2;i++)      // Index pair
		{
			top[i] = i*2;
			bot[i] = i*2+1;
		}
	}

	const int nb1 = CASHSIZE / n;
	assert(CASHSIZE % n == 0);

	#pragma omp parallel
	{
		const int nt = omp_get_num_threads();
		const int id = omp_get_thread_num();
		const int nb = n/nt;
		assert(n % nt == 0);

		int pp, qq;
		bool flg;

		for (int r=0; r<2*nt; r++)
		{
			flg = modulus_pair(2*nt,id,r,pp,qq);

			#pragma omp critical
			cout << "id = " << id << ": r = " << r << ", (pp,qq) = (" << pp << ", " << qq << ")\n";

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
		} // End of r-loop
	} // End of parallel section


	delete [] A;
	delete [] V;
	delete [] top;
	delete [] bot;

	return EXIT_SUCCESS;
}
