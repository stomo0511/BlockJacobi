#ifndef JACOBI_HPP_
#define JACOBI_HPP_

#include <limits>
// Machine epsilon
constexpr double e = std::numeric_limits<double>::epsilon();

void Gen_test_mat(int mode, double cond, const int n, double* D);
void Gen_mat(const int m, const int n, double *A);
void Copy_mat(const int m, const int n, double *A, double *B);
void Show_mat(const int m, const int n, double *A);
void Set_Iden(const int n, double *A);
double Off_d(const int n, double *A);
void music(const int n, int *top, int *bot);
bool modulus_pair(const int w, const int id, const int r, int &p, int &q);

#endif /* JACOBI_HPP_ */
