#ifndef MATRIX_ALGO_H_INCLUDED
#define MATRIX_ALGO_H_INCLUDED

#include "matrix.h"

matrix_type_t matrix_plus(const matrix_type_t a, const matrix_type_t b);
matrix_type_t matrix_minus(const matrix_type_t a, const matrix_type_t b);
matrix_type_t matrix_mul(const matrix_type_t a, const matrix_type_t b);

long double matrix_determinant(const matrix_type_t mat);
matrix_type_t matrix_adjoint(const matrix_type_t mat);
matrix_type_t matrix_inverse(const matrix_type_t mat);

matrix_type_t matrix_transpose(const matrix_type_t mat);

matrix_type_t matrix_diag(unsigned h_w, long double def_val);
matrix_type_t matrix_diag_permut(const matrix_type_t mat);

matrix_type_t matrix_from_two_blocks(const matrix_type_t a, const matrix_type_t b);

long double vector_len(const matrix_type_t vec);

#endif  /* MATRIX_ALGO_H_INCLUDED */
