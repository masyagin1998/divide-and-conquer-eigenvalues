#include "matrix_algo.h"

#include <stddef.h>

matrix_type_t matrix_plus(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j;

    matrix_type_t res = matrix_create(matrix_height(a), matrix_width(a));
    if (res == NULL) {
        goto err0;
    }

    for (i = 1; i <= matrix_height(res); i++) {
        for (j = 1; j <= matrix_width(res); j++) {
            long double cell = matrix_get(a, i, j) + matrix_get(b, i, j);
            matrix_set(res, i, j, cell);
        }
    }

    return res;

 err0:
    return NULL;
}

matrix_type_t matrix_minus(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j;

    matrix_type_t res = matrix_create(matrix_height(a), matrix_width(a));
    if (res == NULL) {
        goto err0;
    }

    for (i = 1; i <= matrix_height(res); i++) {
        for (j = 1; j <= matrix_width(res); j++) {
            long double cell = matrix_get(a, i, j) - matrix_get(b, i, j);
            matrix_set(res, i, j, cell);
        }
    }

    return res;

 err0:
    return NULL;
}

matrix_type_t matrix_mul(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j, k;
    
    matrix_type_t res = matrix_create(matrix_height(a), matrix_width(b));
    if (res == NULL) {
        goto err0;
    }
    
    for (i = 1; i <= matrix_height(res); i++) {
        for (j = 1; j <= matrix_width(res); j++) {
            long double cell = 0.0;
            for (k = 1; k <= matrix_width(a); k++) {
                cell += matrix_get(a, i, k) * matrix_get(b, k, j);
            }
            matrix_set(res, i, j, cell);
        }
    }

    return res;

 err0:
    return NULL;
}
