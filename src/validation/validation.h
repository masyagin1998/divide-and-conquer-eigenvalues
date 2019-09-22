#ifndef VALIDATION_H_INCLUDED
#define VALIDATION_H_INCLUDED

#include "matrix.h"

int matrix_is_square(const matrix_type_t mat);
int matrix_is_symmetric(const matrix_type_t mat, long double eps);
int matrix_is_tridiagonal(const matrix_type_t mat, long double eps);

#endif  /* VALIDATION_H_INCLUDED */
