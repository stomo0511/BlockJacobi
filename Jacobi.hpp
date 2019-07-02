#ifndef JACOBI_HPP_
#define JACOBI_HPP_

// Convergence criterion
#define EPS 1.0E-12

void Gen_mat(const int m, const int n, double *A);
void Copy_mat(const int m, const int n, double *A, double *B);
void Show_mat(const int m, const int n, double *A);
void Set_Iden(const int n, double *A);
double Off_d(const int n, double *A);
void music(const int n, int *top, int *bot);

#endif /* JACOBI_HPP_ */
