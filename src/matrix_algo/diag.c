#include "matrix_algo.h"
#include <stddef.h>

matrix_type_t matrix_diag(unsigned h_w, double def_val)
{
    unsigned i;
    
    matrix_type_t res = matrix_create(h_w, h_w, 0.0);
    if (res == NULL) {
        goto err0;
    }

    for (i = 0; i < h_w; i++) {
        matrix_set(res, i, i, def_val);
    }

    return res;

 err0:
    return NULL;
}

