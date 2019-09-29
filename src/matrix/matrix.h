#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

struct MATRIX;
typedef struct MATRIX* matrix_type_t;

matrix_type_t matrix_create(unsigned h, unsigned w);

void matrix_set_def_val(matrix_type_t mat, long double val);

unsigned matrix_height(const matrix_type_t mat);
unsigned matrix_width(const matrix_type_t mat);

void matrix_set(matrix_type_t mat, unsigned h, unsigned w, long double val);
long double matrix_get(const matrix_type_t mat, unsigned h, unsigned w);

void matrix_free(matrix_type_t mat);

#endif  /* MATRIX_H_INCLUDED */
