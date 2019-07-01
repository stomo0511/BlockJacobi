#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <mkl_cblas.h>

#include "Jacobi.hpp"

using namespace std;

int main(int argc, char **argv)
{
	assert(argc > 2);

	// Matrix size
	const int n = atoi(argv[1]);

	// Submatrix size
	const int nb = atoi(argv[2]);
	assert(n % nb ==0);

	const int w = n / nb;
	assert(w % 2 ==0);

	double *a = new double[n*n];   // Original matrix
	double *oa = new double[n*n];  // Copy of original matrix
	double *v = new double[n*n];   // Right-singular vector matrix
	int *top = new int[w/2];
	int *bot = new int[w/2];

	// Generate matrix
	Gen_mat(n,a);

	#pragma omp parallel
	{
		Copy_mat(n,a,oa);
		Set_Iden(n,v);

		#pragma omp for
		for (int i=0;i<w/2;i++)
		{
			top[i] = i*2;
			bot[i] = i*2+1;
		}
	}

	int k = 0;
	double time = omp_get_wtime();

	double maxt = 1.0;  // convergence criterion
	while (maxt > EPS)
	{
		maxt = 0.0;
		for (int j=0; j<n-1; j++)
		{
			for (int i=0; i<n/2; i++)
			{
				int p = (top[i] > bot[i]) ? top[i] : bot[i];
				int q = (top[i] < bot[i]) ? top[i] : bot[i];

				double *ap = new double[n*nb];
				double *aq = new double[n*nb];

				double x = cblas_ddot(n,a+p*n,1,a+p*n,1);  // x = a_p^T a_p
				double y = cblas_ddot(n,a+q*n,1,a+q*n,1);  // y = a_q^T a_q
				double z = cblas_ddot(n,a+p*n,1,a+q*n,1);  // z = a_p^T a_q

				k++;
			} // End of l-loop
			music(n,top,bot);
		} // End of p-loop
	} // End of while-loop

	time = omp_get_wtime() - time;

	delete [] a;
	delete [] oa;
	delete [] v;
	delete [] top;
	delete [] bot;

	return EXIT_SUCCESS;
}
