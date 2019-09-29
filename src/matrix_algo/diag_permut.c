#include "matrix_algo.h"

#include <assert.h>
#include <stdlib.h>

matrix_type_t matrix_diag_permut(const matrix_type_t mat)
{
    unsigned n, i, j;
    
    matrix_type_t res;
    long double*arr;
    unsigned*idxs;

    assert(matrix_height(mat) == matrix_width(mat));

    n = matrix_height(mat);
    
    arr = (long double*) malloc(n * sizeof(long double));
    if (arr == NULL) {
        goto err0;
    }
    idxs = malloc(n * sizeof(unsigned));
    if (idxs == NULL) {
        goto err1;
    }
    res = matrix_create(n, n);
    if (res == NULL) {
        goto err2;
    }

    for (i = 1; i <= n; i++) {
        arr[i - 1] = matrix_get(mat, i, i);
    }
    for (i = 1; i <= n; i++) {
        idxs[i - 1] = i - 1;
    }
    matrix_set_def_val(res, 0.0);    

    for (i = 1; i <= n; i++) {
        for (j = i + 1; j <= n; j++) {
            if (arr[idxs[i]] > arr[idxs[j]]) {
                unsigned tmp = idxs[i];
                idxs[i] = idxs[j];
                idxs[j] = tmp;
            }
        }
    }

    for (i = 1; i <= matrix_height(res); i++) {
        matrix_set(res, i, idxs[i - 1], 1.0);
    }

    free(arr);
    free(idxs);

    return res;

 err2:
    free(idxs);
 err1:
    free(arr);
 err0:
    return NULL;
}
