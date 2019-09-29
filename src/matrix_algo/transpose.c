#include "matrix_algo.h"

#include <stddef.h>

matrix_type_t matrix_transpose(const matrix_type_t mat)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(matrix_width(mat), matrix_height(mat));
    if (res == NULL) {
        goto err0;
    }

    for (i = 1; i <= matrix_width(mat); i++) {
        for (j = 1; j <= matrix_height(mat); j++) {
            long double cell = matrix_get(mat, j, i);
            matrix_set(res, i, j, cell);
        }
    }

    return res;

 err0:
    return NULL;
}
