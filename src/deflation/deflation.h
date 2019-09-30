#ifndef DEFLATION_H_INCLUDED
#define DEFLATION_H_INCLUDED

#include "matrix.h"

int deflate(const matrix_type_t D, matrix_type_t*v, matrix_type_t*v_prime, matrix_type_t*eigenvalues, matrix_type_t*eigenvectors, unsigned*n_deflated, matrix_type_t*G, long double eps);

#endif  /* DEFLATION_H_INCLUDED */
