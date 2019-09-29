#include "matrix_algo.h"

#include <math.h>

long double vector_norm(const matrix_type_t vec)
{
    unsigned i;
    long double len = 0.0;
    
    if (matrix_height(vec) == 1) {
        for (i = 1; i <= matrix_width(vec); i++) {
            len += (matrix_get(vec, 1, i) * matrix_get(vec, 1, i));
        }
    } else if (matrix_width(vec) == 1) {
        for (i = 1; i <= matrix_height(vec); i++) {
            len += (matrix_get(vec, i, 1) * matrix_get(vec, i, 1));
        }
    }

    return sqrtl(len);
}
