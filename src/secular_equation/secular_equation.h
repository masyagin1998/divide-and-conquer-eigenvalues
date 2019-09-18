#ifndef SECULAR_EQUATION_H_INCLUDED
#define SECULAR_EQUATION_H_INCLUDED

#include "matrix.h"

double solve_secular_equation_d0_plus_inf(const matrix_type_t D, const matrix_type_t u, double ro, double eps);
double solve_secular_equation_common(const matrix_type_t D, const matrix_type_t u, double ro, unsigned i, double eps);
double solve_secular_equation_minus_inf_dn(const matrix_type_t D, const matrix_type_t u, double ro, double eps);

#endif  /* SECULAR_EQUATION_H_INCLUDED */
