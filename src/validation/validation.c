#include "validation.h"

#include <math.h>

__inline__ int matrix_is_square(const matrix_type_t mat)
{
    return (matrix_height(mat) == matrix_width(mat));
}

#define are_equal(v1, v2) ((fabs((v1) - (v2))) < eps)

int matrix_is_symmetric(const matrix_type_t mat, double eps)
{
    unsigned i, j;
    
    if (!matrix_is_square(mat)) {
        return 0;
    }

    for (i = 0; i < matrix_height(mat); i++) {
        for (j = 0; j < (i + 1); j++) {
            if (!are_equal(matrix_get(mat, i, j), matrix_get(mat, j, i))) {
                return 0;
            }
        }
    }

    return 1;
}

#define is_zero(v) ((v) < eps)

int matrix_is_tridiagonal(const matrix_type_t mat, double eps)
{
    unsigned i, j;
    
    if (!matrix_is_square(mat)) {
        return 0;
    }

    for (i = 0; i < matrix_height(mat); i++) {
        for (j = 0; j < matrix_width(mat); j++) {
            double cell = matrix_get(mat, i, j);
            if ((i == j) || (i - 1 == j) || (i + 1 == j)) {
                if (is_zero(cell)) {
                    return 0;
                }
            } else {
                if (!is_zero(cell)) {
                    return 0;
                }
            }
        }
    }

    return 1;
}
