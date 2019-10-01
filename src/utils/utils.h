#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include "matrix.h"

/*
  Format of matrix file:
  two positivie values in first string: height and width;
  height string with width long double values in every string;
 */

void matrix_print(const matrix_type_t mat, long double eps);
matrix_type_t read_matrix_from_file(const char*fname);
matrix_type_t matrix_from_array(unsigned h, unsigned w, const long double arr[h][w]);
int save_matrix_to_file(const matrix_type_t mat, const char*fname);
matrix_type_t matrix_copy(const matrix_type_t mat);

#endif  /* UTILS_H_INCLUDED */
