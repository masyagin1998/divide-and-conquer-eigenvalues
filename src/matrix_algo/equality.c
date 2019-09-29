#include "matrix_algo.h"
#include "macro.h"

#include <stddef.h>

int matrices_are_equal(const matrix_type_t a, const matrix_type_t b, long double eps)
{
    unsigned i, j;
    
    if (matrix_height(a) != matrix_height(b)) {
        return 0;
    }

    if (matrix_width(a) != matrix_width(b)) {
        return 0;
    }

    for (i = 1; i <= matrix_height(a); i++) {
        for (j = 1; j <= matrix_width(a); j++) {
            if (!are_equal(matrix_get(a, i, j), matrix_get(b, i, j))) {
                return 0;
            }
        }
    }

    return 1;
}
