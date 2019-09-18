#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include "matrix.h"

/*
  Format of matrix file:
  two positivie values in first string: height and width;
  height string with width double values in every string;
 */

void matrix_print_part(const matrix_type_t mat, unsigned h_min, unsigned w_min, unsigned h_max, unsigned w_max);
void matrix_print(const matrix_type_t mat);

matrix_type_t read_matrix_from_file(const char*fname);
int save_matrix_to_file(const matrix_type_t, const char*fname);

#endif  /* UTILS_H_INCLUDED */
