/*
 * Jacobi.hpp
 *
 *  Created on: 2019/06/12
 *      Author: stomo
 */

#ifndef JACOBI_HPP_
#define JACOBI_HPP_

// Convergence criterion
#define EPS 1.0E-12

void Gen_symmat(const int n, double *a);
void Gen_mat(const int n, double *a);
void Copy_symmat(const int n, double *a, double *b);
void Copy_mat(const int n, double *a, double *b);
void Set_Iden(const int n, double *a);
double Off_d(const int n, double *a);
void Sym_schur2(const int n, double *a, const int p, const int q, double *c, double *s);
void Search_max(const int n, double *a, int *p, int *q);
void PSearch_max(const int n, double *a, int *p, int *q);
void Givens(const int n, double *a, const int p, const int q, const double c, const double s);
void GivensR(const int n, double *v, const int p, const int q, const double c, const double s);
void GivensL(const int n, double *v, const int p, const int q, const double c, const double s);
void Givens2(const int n, double *a, const int p, const int q, const double c, const double s, double *b);
void GivensR2(const int n, double *v, const int p, const int q, const double c, const double s, double *u);
void GivensL2(const int n, double *v, const int p, const int q, const double c, const double s, double *u);
void music(const int n, int *top, int *bot);
double Residure(const int n, double* oa, double* a, double *v);

#endif /* JACOBI_HPP_ */
