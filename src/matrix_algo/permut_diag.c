#include "matrix_algo.h"

#include <stdlib.h>

matrix_type_t matrix_diag_permut(const matrix_type_t mat)
{
    unsigned i, j;
    
    matrix_type_t res;
    double*arr;
    unsigned*idxs;
    
    arr = (double*) malloc(matrix_height(mat) * sizeof(double));
    if (arr == NULL) {
        goto err0;
    }
    for (i = 0; i < matrix_height(mat); i++) {
        arr[i] = matrix_get(mat, i, i);
    }
    idxs = malloc(matrix_height(mat) * sizeof(unsigned));
    if (idxs == NULL) {
        goto err1;
    }
    for (i = 0; i < matrix_height(mat); i++) {
        idxs[i] = i;
    }

    for (i = 0; i < matrix_height(mat); i++) {
        for (j = i + 1; j < matrix_height(mat); j++) {
            if (arr[idxs[i]] > arr[idxs[j]]) {
                unsigned tmp = idxs[i];
                idxs[i] = idxs[j];
                idxs[j] = tmp;
            }
        }
    }

    for (i = 0; i < matrix_height(mat) / 2; i++) {
        unsigned tmp = idxs[i];
        idxs[i] = idxs[matrix_height(mat) - 1 - i];
        idxs[matrix_height(mat) - 1 - i] = tmp;
    }

    res = matrix_create(matrix_height(mat), matrix_width(mat), 0.0);
    if (res == NULL) {
        goto err2;
    }
    for (i = 0; i < matrix_height(res); i++) {
        matrix_set(res, i, idxs[i], 1.0);
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
