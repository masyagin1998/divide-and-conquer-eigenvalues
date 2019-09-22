#include "matrix_algo.h"

#include <stddef.h>

matrix_type_t matrix_from_two_blocks(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(matrix_height(a) + matrix_height(b), matrix_width(a) + matrix_width(b), 0.0);
    if (res == NULL) {
        goto err0;
    }
    
    for (i = 0; i < matrix_height(a); i++) {
        for (j = 0; j < matrix_width(a); j++) {
            matrix_set(res, i, j, matrix_get(a, i, j));
        }
    }

    for (i = 0; i < matrix_height(b); i++) {
        for (j = 0; j < matrix_width(b); j++) {
            matrix_set(res, i + matrix_height(a), j + matrix_width(a), matrix_get(b, i, j));
        }
    }

    return res;

 err0:
    return NULL;
}
