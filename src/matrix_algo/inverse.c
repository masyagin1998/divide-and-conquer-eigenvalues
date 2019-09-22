#include "matrix_algo.h"

#include <stddef.h>

void matrix_cofactor(const matrix_type_t mat, matrix_type_t tmp, unsigned p, unsigned q)
{
    unsigned i = 0, j = 0;
    unsigned row, col;

    for (row = 0; row < matrix_height(mat); row++) {
        for (col = 0; col < matrix_width(mat); col++) {
            if ((row != p) && (col != q)) {
                matrix_set(tmp, i, j, matrix_get(mat, row, col));
                j++;
                if (j == matrix_height(mat) - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

long double matrix_determinant(const matrix_type_t mat)
{
    unsigned i;
    
    long double det;
    matrix_type_t tmp;
    int sign;

    if (matrix_height(mat) == 1) {
        return matrix_get(mat, 0, 0);
    }

    det = 0.0;
    tmp = matrix_create(matrix_height(mat) - 1, matrix_width(mat) - 1, 0.0);
    sign = 1;

    for (i = 0; i < matrix_height(mat); i++) {
        matrix_cofactor(mat, tmp, 0, i);
        det += sign * matrix_get(mat, 0, i) * matrix_determinant(tmp);
        sign = -sign;
    }

    matrix_free(tmp);

    return det;
}

matrix_type_t matrix_adjoint(const matrix_type_t mat)
{
    unsigned i, j;
    matrix_type_t res;
    matrix_type_t tmp;
    int sign;

    res = matrix_create(matrix_height(mat), matrix_width(mat), 0.0);
    if (res == NULL) {
        goto err0;
    }

    if (matrix_height(res) == 1) {
        matrix_set(res, 0, 0, matrix_get(mat, 0, 0));
        return res;
    }

    tmp = matrix_create(matrix_height(mat) - 1, matrix_width(mat) - 1, 0.0);
    if (tmp == NULL) {
        goto err1;
    }

    sign = 1;

    for (i = 0; i < matrix_height(mat); i++) {
        for (j = 0; j < matrix_width(mat); j++) {
            long double cell;
            
            matrix_cofactor(mat, tmp, i, j);

            sign = ((i+j)%2==0)? 1: -1;

            cell = (sign) * matrix_determinant(tmp);
            matrix_set(res, j, i, cell);
        }
    }

    matrix_free(tmp);

    return res;

 err1:
    matrix_free(res);
 err0:
    return NULL;
}

matrix_type_t matrix_inverse(const matrix_type_t mat)
{
    unsigned i, j;
    matrix_type_t res;
    matrix_type_t adj;
    long double det;

    res = matrix_create(matrix_height(mat), matrix_width(mat), 0.0);
    if (res == NULL) {
        goto err0;
    }

    adj = matrix_adjoint(mat);
    if (adj == NULL) {
        goto err1;
    }

    det = matrix_determinant(mat);

    for (i = 0; i < matrix_height(res); i++) {
        for (j = 0; j < matrix_width(res); j++) {
            long double cell = matrix_get(adj, i, j) / det;
            matrix_set(res, i, j, cell);
        }
    }

    matrix_free(adj);

    return res;

 err1:
    matrix_free(res);
 err0:
    return NULL;
}
