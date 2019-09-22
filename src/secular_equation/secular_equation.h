#ifndef SECULAR_EQUATION_H_INCLUDED
#define SECULAR_EQUATION_H_INCLUDED

#include "matrix.h"

long double solve_secular_equation_d0_plus_inf(const matrix_type_t D, const matrix_type_t u, long double ro, long double eps);
long double solve_secular_equation_common(const matrix_type_t D, const matrix_type_t u, long double ro, unsigned i, long double eps);
long double solve_secular_equation_minus_inf_dn(const matrix_type_t D, const matrix_type_t u, long double ro, long double eps);

#endif  /* SECULAR_EQUATION_H_INCLUDED */
