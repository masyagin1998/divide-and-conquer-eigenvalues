#include <stdio.h>

#include "utils.h"
#include "validation.h"
#include "divide_and_conquer.h"

#define EPS 0.0000001

int main(int argc, char**argv)
{
    const char*error_str = NULL;
    int        error_code = 0;
    
    matrix_type_t mat;

    if (argc != 2) {
        printf("usage:   %s path/to/input/matrix\n", argv[0]);
        printf("example: %s data/in.mat\n", argv[0]);
        return -1;
    }

    mat = read_matrix_from_file(argv[1]);
    if (mat == NULL) {
        error_code = -2;
        error_str = "unable to read matrix from file";
        goto err0;
    }

    if (!matrix_is_square(mat)) {
        error_code = -3;
        error_str = "matrix is not square";
        goto err1;
    }

    if (!matrix_is_symmetric(mat, EPS)) {
        error_code = -4;
        error_str = "matrix is not symmetric";
        goto err1;
    }

    if (!matrix_is_tridiagonal(mat, EPS)) {
        error_code = -5;
        error_str = "matrix is not tridiagonal";
        goto err1;
    }

    printf("INPUT MATRIX:\n");
    matrix_print(mat);
    matrix_divide_and_conquer(mat, EPS);

    matrix_free(mat);

    return 0;
 err1:
    matrix_free(mat);
 err0:
    printf("error: %s\n", error_str);
    return error_code;
}
