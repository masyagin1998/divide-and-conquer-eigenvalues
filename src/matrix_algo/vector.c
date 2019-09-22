#include "matrix_algo.h"

#include <math.h>

long double vector_len(const matrix_type_t vec)
{
    unsigned i;
    long double len = 0.0;
    
    if (matrix_height(vec) == 1) {
        for (i = 0; i < matrix_width(vec); i++) {
            len += (matrix_get(vec, 0, i) * matrix_get(vec, 0, i));
        }
    } else if (matrix_width(vec) == 1) {
        for (i = 0; i < matrix_height(vec); i++) {
            len += (matrix_get(vec, i, 0) * matrix_get(vec, i, 0));
        }
    }

    return sqrt(len);
}
