#include "matrix_algo.h"

#include <assert.h>
#include <stddef.h>

void matrix_cofactor(const matrix_type_t mat, matrix_type_t tmp, unsigned p, unsigned q)
{
    unsigned n = matrix_height(mat), i = 1, j = 1;
    unsigned row, col;

    assert(matrix_height(mat) == matrix_width(mat));
    assert(matrix_height(tmp) == matrix_width(tmp));

    for (row = 1; row <= n; row++) {
        for (col = 1; col <= n; col++) {
            if ((row != p) && (col != q)) {
                matrix_set(tmp, i, j, matrix_get(mat, row, col));
                j++;
                if (j == n) {
                    j = 1;
                    i++;
                }
            }
        }
    }
}

int matrix_determinant(const matrix_type_t mat, long double*det)
{
    int r;
    unsigned n, i;
    
    matrix_type_t tmp;
    int sign;

    assert(matrix_height(mat) == matrix_width(mat));

    n = matrix_height(mat);

    if (n == 1) {
        return matrix_get(mat, 0, 0);
    }

    tmp = matrix_create(n - 1, n - 1);
    if (tmp == NULL) {
        r = -1;
        goto err0;
    }

    (*det) = 0.0;
    sign = 1;

    for (i = 1; i <= n; i++) {
        long double det1;
        r = matrix_determinant(tmp, &det1);
        if (r) {
            goto err1;
        }
        matrix_cofactor(mat, tmp, 1, i);
        (*det) += sign * matrix_get(mat, 1, i) * det1;
        sign = -sign;
    }

    matrix_free(tmp);

 err1:
    matrix_free(tmp);
 err0:
    return r;
}

matrix_type_t matrix_adjoint(const matrix_type_t mat)
{
    int r;
    
    unsigned i, j, n;
    matrix_type_t res;
    matrix_type_t tmp;
    int sign;

    assert(matrix_height(mat) == matrix_width(mat));

    n = matrix_height(mat);

    res = matrix_create(n, n);
    if (res == NULL) {
        goto err0;
    }

    if (n == 1) {
        long double cell = matrix_get(mat, 1, 1);
        matrix_set(res, 1, 1, cell);
        return res;
    }

    tmp = matrix_create(matrix_height(mat) - 1, matrix_width(mat) - 1);
    if (tmp == NULL) {
        goto err1;
    }    

    sign = 1;

    for (i = 1; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            long double det;
            long double cell;
            
            matrix_cofactor(mat, tmp, i, j);
            r = matrix_determinant(tmp, &det);
            if (r) {
                goto err2;
            }

            sign = ((i + j) % 2 == 0) ? 1: -1;

            cell = sign * det;
            matrix_set(res, j, i, cell);
        }
    }

    matrix_free(tmp);

    return res;

 err2:
    matrix_free(tmp);
 err1:
    matrix_free(res);
 err0:
    return NULL;
}

matrix_type_t matrix_inverse(const matrix_type_t mat)
{
    int r;
    unsigned i, j;
    matrix_type_t res;
    matrix_type_t adj;
    long double det;

    res = matrix_create(matrix_height(mat), matrix_width(mat));
    if (res == NULL) {
        goto err0;
    }
    adj = matrix_adjoint(mat);
    if (adj == NULL) {
        goto err1;
    }
    r = matrix_determinant(mat, &det);
    if (r) {
        goto err2;
    }

    for (i = 0; i < matrix_height(res); i++) {
        for (j = 0; j < matrix_width(res); j++) {
            long double cell = matrix_get(adj, i, j) / det;
            matrix_set(res, i, j, cell + 0.0);
        }
    }

    matrix_free(adj);

    return res;

 err2:
    matrix_free(adj);
 err1:
    matrix_free(res);
 err0:
    return NULL;
}
