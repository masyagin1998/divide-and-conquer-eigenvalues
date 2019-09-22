#include "matrix_algo.h"

#include <stddef.h>

matrix_type_t matrix_transpose(const matrix_type_t mat)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(matrix_width(mat), matrix_height(mat), 0.0);
    if (res == NULL) {
        goto err0;
    }

    for (i = 0; i < matrix_width(mat); i++) {
        for (j = 0; j < matrix_height(mat); j++) {
            matrix_set(res, i, j, matrix_get(mat, j, i));
        }
    }

    return res;

 err0:
    return NULL;
}
