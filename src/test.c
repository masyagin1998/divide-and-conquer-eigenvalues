/*
  If TEST if defined, program will test basic functionality
  from src/matrix, src/macro, src/matrix_algo, src/utils,
  else program will run divide-and-conquer eigenvalues algorithm.
*/
/* #define TEST */

#ifdef TEST

#include "matrix.h"
#include "utils.h"
#include "validation.h"
#include "matrix_algo.h"

#include <assert.h>

#define EPS 0.000001

int main()
{
    {
        /* matrix_plus */
        const long double a_arr[3][3] =
        {
         { 1.0, 2.0, 3.0 },
         { 4.0, 5.0, 6.0 },
         { 7.0, 8.0, 9.0 },
        };
        matrix_type_t a = matrix_from_array(3, 3, a_arr);
        const long double b_arr[3][3] =
        {
         { 9.0, 8.0, 7.0 },
         { 6.0, 5.0, 4.0 },
         { 3.0, 2.0, 1.0 },
        };
        matrix_type_t b = matrix_from_array(3, 3, b_arr);
        const long double res_arr[3][3] =
        {
         { 10.0, 10.0, 10.0 },
         { 10.0, 10.0, 10.0 },
         { 10.0, 10.0, 10.0 },
        };
        matrix_type_t true_res = matrix_from_array(3, 3, res_arr);

        matrix_type_t res = matrix_plus(a, b);

        assert(matrices_are_equal(res, true_res, EPS));
        matrix_free(res);
    }

    {
        /* matrix_minus */
        const long double a_arr[3][3] =
        {
         { 1.0, 2.0, 3.0 },
         { 4.0, 5.0, 6.0 },
         { 7.0, 8.0, 9.0 },
        };
        matrix_type_t a = matrix_from_array(3, 3, a_arr);
        const long double b_arr[3][3] =
        {
         { 9.0, 8.0, 7.0 },
         { 6.0, 5.0, 4.0 },
         { 3.0, 2.0, 1.0 },
        };
        matrix_type_t b = matrix_from_array(3, 3, b_arr);
        const long double res_arr[3][3] =
        {
         { -8.0, -6.0, -4.0 },
         { -2.0,  0.0,  2.0 },
         {  4.0,  6.0,  8.0 },
        };
        matrix_type_t true_res = matrix_from_array(3, 3, res_arr);

        matrix_type_t res = matrix_minus(a, b);

        assert(matrices_are_equal(res, true_res, EPS));
        matrix_free(res);
    }

    {
        /* matrix_mul */
    }

    {
        /* */
    }

    {
        /* */
    }

    {
        /* */
    }

    {
        /* */
    }    
}

#else

#include "utils.h"
#include "validation.h"
#include "matrix_algo.h"
#include "divide_and_conquer.h"

#include <stdio.h>

#define EPS 0.00001

int main(int argc, char**argv)
{
    int r;
    const char*error_str = NULL;
    int        error_code = 0;
    
    matrix_type_t T, Q, L, Q_t, tmp1, tmp2;

    if (argc != 2) {
        printf("usage:   %s path/to/input/matrix\n", argv[0]);
        printf("example: %s data/in.mat\n", argv[0]);
        return -1;
    }

    T = read_matrix_from_file(argv[1]);
    if (T == NULL) {
        error_code = -2;
        error_str = "unable to read matrix from file";
        goto err0;
    }

    if (!matrix_is_square(T)) {
        error_code = -3;
        error_str = "matrix is not square";
        goto err1;
    }

    if (!matrix_is_symmetric(T, EPS)) {
        error_code = -4;
        error_str = "matrix is not symmetric";
        goto err1;
    }

    if (!matrix_is_tridiagonal(T, EPS)) {
        error_code = -5;
        error_str = "matrix is not tridiagonal";
        goto err1;
    }

    r = matrix_divide_and_conquer(T, &Q, &L, EPS);
    if (r) {
        error_code = r;
        error_str = "unable to run divide-and-conquer algorithm";
        goto err1;
    }
    printf("\n");
    printf("in:\n");
    matrix_print(T);
    printf("Q:\n");
    matrix_print(Q);
    printf("L:\n");
    matrix_print(L);
    Q_t = matrix_transpose(Q);
    printf("Q t:\n");
    matrix_print(Q_t);
    printf("Q * L * Q t:\n");
    tmp1 = matrix_mul(Q, L);
    tmp2 = matrix_mul(tmp1, Q_t);
    matrix_print(tmp2);
    matrix_free(tmp1);
    matrix_free(tmp2);
    matrix_free(Q_t);

    return 0;
 err1:
    matrix_free(T);
 err0:
    printf("error: %s\n", error_str);
    return error_code;
}

#endif  /* TEST */
