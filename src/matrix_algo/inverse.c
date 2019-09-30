#include "matrix_algo.h"

#include <stddef.h>

void matrix_cofactor(const matrix_type_t mat, matrix_type_t tmp, unsigned p, unsigned q)
{
    unsigned i = 1, j = 1;
    unsigned row, col;

    for (row = 1; row <= matrix_height(mat); row++) {
        for (col = 1; col <= matrix_width(mat); col++) {
            if ((row != p) && (col != q)) {
                matrix_set(tmp, i, j, matrix_get(mat, row, col));
                j++;
                if (j == matrix_height(mat)) {
                    j = 1;
                    i++;
                }
            }
        }
    }
}

double matrix_determinant(const matrix_type_t mat)
{
    unsigned i;
    
    double det;
    matrix_type_t tmp;
    int sign;

    if (matrix_height(mat) == 1) {
        return matrix_get(mat, 1, 1);
    }

    det = 0.0;
    tmp = matrix_create(matrix_height(mat) - 1, matrix_width(mat) - 1);
    matrix_set_def_val(tmp, 0.0);
    sign = 1;

    for (i = 1; i <= matrix_height(mat); i++) {
        matrix_cofactor(mat, tmp, 1, i);
        det += sign * matrix_get(mat, 1, i) * matrix_determinant(tmp);
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

    res = matrix_create(matrix_height(mat), matrix_width(mat));
    if (res == NULL) {
        goto err0;
    }
    matrix_set_def_val(res, 0.0);
    if (matrix_height(res) == 1) {
        matrix_set(res, 1, 1, 1.0);
        return res;
    }

    tmp = matrix_create(matrix_height(mat) - 1, matrix_width(mat) - 1);    
    if (tmp == NULL) {
        goto err1;
    }
    matrix_set_def_val(tmp, 0.0);

    sign = 1;

    for (i = 1; i <= matrix_height(mat); i++) {
        for (j = 1; j <= matrix_width(mat); j++) {
            double cell;
            
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
    double det;

    res = matrix_create(matrix_height(mat), matrix_width(mat));
    if (res == NULL) {
        goto err0;
    }
    matrix_set_def_val(res, 0.0);

    adj = matrix_adjoint(mat);
    if (adj == NULL) {
        goto err1;
    }

    det = matrix_determinant(mat);

    for (i = 1; i <= matrix_height(res); i++) {
        for (j = 1; j <= matrix_width(res); j++) {
            double cell = matrix_get(adj, i, j) / det;
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
