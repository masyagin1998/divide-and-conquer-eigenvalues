#include "matrix.h"

#include <assert.h>
#include <stdlib.h>

struct MATRIX
{
    unsigned h;    
    unsigned w;    
    long double vals[];
};

matrix_type_t matrix_create(unsigned h, unsigned w)
{
    matrix_type_t mat;

    assert(h >= 1);
    assert(w >= 1);    
    
    mat = (matrix_type_t) malloc(sizeof(struct MATRIX) + ((h * w) * sizeof(long double)));
    if (mat == NULL) {
        goto err0;
    }
    
    mat->h = h;
    mat->w = w;

    return mat;

 err0:
    return NULL;
}

void matrix_set_def_val(matrix_type_t mat, long double val)
{
    unsigned i;
    for (i = 0; i < (mat->h * mat->w); i++) {
        mat->vals[i] = val;
    }
}

__inline__ unsigned matrix_height(const matrix_type_t mat)
{
    return mat->h;
}

__inline__ unsigned  matrix_width(const matrix_type_t mat)
{
    return mat->w;    
}

__inline__ void matrix_set(matrix_type_t mat, unsigned h, unsigned w, long double val)
{
    assert((h >= 1) && (h <= mat->h));
    assert((w >= 1) && (w <= mat->w));

    mat->vals[((h - 1) * mat->w) + (w - 1)] = val;
}

__inline__ long double matrix_get(const matrix_type_t mat, unsigned h, unsigned w)
{
    assert((h >= 1) && (h <= mat->h));
    assert((w >= 1) && (w <= mat->w));

    return mat->vals[((h - 1) * mat->w) + (w - 1)];
}

void matrix_free(matrix_type_t mat)
{
    free(mat);
}
