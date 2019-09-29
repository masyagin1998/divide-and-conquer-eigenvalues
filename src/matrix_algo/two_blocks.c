#include "matrix_algo.h"

#include <stddef.h>

matrix_type_t matrix_from_two_blocks(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(matrix_height(a) + matrix_height(b), matrix_width(a) + matrix_width(b));
    if (res == NULL) {
        goto err0;
    }

    matrix_set_def_val(res, 0.0);
    
    for (i = 1; i <= matrix_height(a); i++) {
        for (j = 1; j <= matrix_width(a); j++) {
            matrix_set(res, i, j, matrix_get(a, i, j));
        }
    }

    for (i = 1; i <= matrix_height(b); i++) {
        for (j = 1; j <= matrix_width(b); j++) {
            matrix_set(res, i + matrix_height(a), j + matrix_width(a), matrix_get(b, i, j));
        }
    }

    return res;

 err0:
    return NULL;
}
