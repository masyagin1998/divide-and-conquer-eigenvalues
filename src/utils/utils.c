#include "matrix.h"

#include <stdio.h>
#include <math.h>

#define MAT_VAL "%5.5Lf"

#define MAT_VAL_LEN 7

void matrix_print(const matrix_type_t mat, long double eps)
{
    unsigned i, j;
    long double cell;

    for (i = 1; i <= matrix_height(mat); i++) {
        printf("|");
        for (j = 1; j <= matrix_width(mat) - 1; j++) {
            cell = matrix_get(mat, i, j);
            if (fabsl(cell) <= eps) {
                cell = 0.0;
            }
            printf(MAT_VAL" ", cell);
        }
        cell = matrix_get(mat, i, matrix_width(mat));
        if (fabsl(cell) <= eps) {
            cell = 0.0;
        }    
        printf(MAT_VAL"|\n", cell);
    }
}

matrix_type_t read_matrix_from_file(const char*fname)
{
    matrix_type_t mat;
    FILE*f;
    unsigned h, w;
    unsigned i, j;
    
    f = fopen(fname, "rb");
    if (f == NULL) {
        goto err0;
    }

    if (fscanf(f, "%u %u\n", &h, &w) != 2) {
        goto err1;
    }

    mat = matrix_create(h, w);
    if (mat == NULL) {
        goto err1;
    }

    for (i = 1; i <= h; i++) {
        for (j = 1; j <= w; j++) {
            long double cell;
            if (fscanf(f, "%Lf", &cell) != 1) {
                goto err2;
            }
            matrix_set(mat, i, j, cell);
        }
    }

    fclose(f);

    return mat;

 err2:
    matrix_free(mat);
 err1:
    fclose(f);
    /* remove(fname); */
 err0:
    return NULL;
}

matrix_type_t matrix_from_array(unsigned h, unsigned w, const long double arr[h][w])
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(h, w);
    if (res == NULL) {
        goto err0;
    }

    for (i = 1; i <= h; i++) {
        for (j = 1; j <= w; j++) {
            long double cell = arr[i - 1][j - 1];
            matrix_set(res, i, j, cell);
        }
    }

    return res;

 err0:
    return NULL;
}

int save_matrix_to_file(const matrix_type_t mat, const char*fname)
{
    int r;
    FILE*f;
    unsigned i, j;

    f = fopen(fname, "wb");
    if (f == NULL) {
        r = -1;
        goto err0;
    }

    fprintf(f, "%u %u\n", matrix_height(mat), matrix_width(mat));

    for (i = 1; i <= matrix_height(mat); i++) {
        for (j = 1; j <= matrix_width(mat) - 1; j++) {
            fprintf(f, "%Lf ", matrix_get(mat, i, j));
        }
        fprintf(f, "%Lf\n", matrix_get(mat, i, matrix_width(mat) - 1));
    }

    return 0;

 err0:
    return r;
}

matrix_type_t matrix_copy(const matrix_type_t mat)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(matrix_height(mat), matrix_width(mat));
    if (res == NULL) {
        goto err0;
    }

    for (i = 1; i <= matrix_height(mat); i++) {
        for (j = 1; j <= matrix_width(mat); j++) {
            long double cell = matrix_get(mat, i, j);
            matrix_set(res, i, j, cell);
        }
    }

    return res;

 err0:
    return NULL;
}
