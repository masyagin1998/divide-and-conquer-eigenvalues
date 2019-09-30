#ifndef SECULAR_EQUATION_H_INCLUDED
#define SECULAR_EQUATION_H_INCLUDED

#include "matrix.h"

long double solve_secular_equation(long double rho, const matrix_type_t D, const matrix_type_t v_prime, unsigned i, long double lambda_init, unsigned n_deflated, long double eps);

#endif  /* SECULAR_EQUATION_H_INCLUDED */
