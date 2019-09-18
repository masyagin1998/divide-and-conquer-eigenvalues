#include "matrix.h"

#include <stdlib.h>

struct MATRIX
{
    double**vals;
    unsigned w;    
    unsigned h;
};

matrix_type_t matrix_create(unsigned h, unsigned w, double def_val)
{
    unsigned i, j;
    
    matrix_type_t mat = (matrix_type_t) malloc(sizeof(struct MATRIX));
    if (mat == NULL) {
        goto err0;
    }
    mat->h = h;
    mat->w = w;

    mat->vals = (double**) malloc(mat->h * sizeof(double*));
    if (mat->vals == NULL) {
        goto err1;
    }
    for (i = 0; i < mat->h; i++) {
        mat->vals[i] = (double*) malloc(mat->w * sizeof(double));
        if (mat->vals[i] == NULL) {
            goto err2;
        }
    }

    for (i = 0; i < mat->h; i++) {
        for (j = 0; j < mat->w; j++) {
            mat->vals[i][j] = def_val;
        }
    }

    return mat;

 err2:
    for (j = 0; j < i; j++) {
        free(mat->vals[j]);
    }
 err1:
    free(mat);
 err0:
    return NULL;
}

__inline__ unsigned matrix_height(const matrix_type_t mat)
{
    return mat->h;
}

__inline__ unsigned  matrix_width(const matrix_type_t mat)
{
    return mat->w;    
}

__inline__ void matrix_set(matrix_type_t mat, unsigned h, unsigned w, double val)
{
    mat->vals[h][w] = val;
}

__inline__ double matrix_get(const matrix_type_t mat, unsigned h, unsigned w)
{
    return mat->vals[h][w];
}

void matrix_free(matrix_type_t mat)
{
    unsigned i;
    for (i = 0; i < mat->w; i++) {
        free(mat->vals[i]);
    }
    free(mat->vals);
    free(mat);
}
