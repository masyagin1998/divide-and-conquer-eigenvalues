#include "divide_and_conquer.h"
#include "secular_equation.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>

static matrix_type_t matrix_from_two_blocks(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(matrix_height(a) + matrix_height(b), matrix_width(a) + matrix_width(b), 0.0);

    for (i = 0; i < matrix_height(a); i++) {
        for (j = 0; j < matrix_width(a); j++) {
            matrix_set(res, i, j, matrix_get(a, i, j));
        }
    }

    for (i = 0; i < matrix_height(b); i++) {
        for (j = 0; j < matrix_width(b); j++) {
            matrix_set(res, i + matrix_height(a), j + matrix_width(a), matrix_get(b, i, j));
        }
    }

    return res;
}

static matrix_type_t matrix_transpose(const matrix_type_t mat)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(matrix_width(mat), matrix_height(mat), 0.0);

    for (i = 0; i < matrix_width(mat); i++) {
        for (j = 0; j < matrix_height(mat); j++) {
            matrix_set(res, i, j, matrix_get(mat, j, i));
        }
    }

    return res;
}

static matrix_type_t matrix_mul(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j, k;
    
    matrix_type_t res = matrix_create(matrix_height(a), matrix_width(b), 0.0);
    for (i = 0; i < matrix_height(res); i++) {
        for (j = 0; j < matrix_width(res); j++) {
            double cell = 0.0;
            for (k = 0; k < matrix_width(a); k++) {
                cell += matrix_get(a, i, k) * matrix_get(b, k, j);
            }
            matrix_set(res, i, j, cell);
        }
    }
    

    return res;
}

#include <stdio.h>

static matrix_type_t matrix_get_permut(const matrix_type_t mat)
{
    unsigned i, j;
    
    matrix_type_t res;
    double*arr;
    unsigned*idxs;
    
    arr = (double*) malloc(matrix_height(mat) * sizeof(double));
    for (i = 0; i < matrix_height(mat); i++) {
        arr[i] = matrix_get(mat, i, i);
    }
    idxs = malloc(matrix_height(mat) * sizeof(unsigned));
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
    for (i = 0; i < matrix_height(res); i++) {
        matrix_set(res, i, idxs[i], 1.0);
    }

    return res;
}

static void matrix_divide_and_conquer_inner(matrix_type_t mat, unsigned b_ind, unsigned e_ind, matrix_type_t*Q, matrix_type_t*lambda, double eps)
{
    unsigned m;
    double b_m;
    unsigned i;

    matrix_type_t tmp;

    matrix_type_t Q1;
    matrix_type_t lambda1;

    matrix_type_t Q2;
    matrix_type_t lambda2;

    matrix_type_t D;
    matrix_type_t P;
    matrix_type_t P_t;

    matrix_type_t Q_t;
    matrix_type_t v;
    matrix_type_t u;

    double*lambdas;
    
    if ((e_ind - b_ind == 1)) {
        (*lambda) = matrix_create(1, 1, matrix_get(mat, b_ind, b_ind));
        (*Q) = matrix_create(1, 1, 1.0);
        return;
    }

    /* m and b_m. */
    m = ((e_ind - b_ind) / 2) + b_ind;
    b_m = matrix_get(mat, m, m - 1);

    /* Create T1. */
    matrix_set(mat, m - 1, m - 1, matrix_get(mat, m - 1, m - 1) - b_m);
    /* Create T2. */
    matrix_set(mat, m, m, matrix_get(mat, m, m) - b_m);

    /* Process T1 and get Q1 and λ1. */
    matrix_divide_and_conquer_inner(mat, b_ind, m, &Q1, &lambda1, eps);
    /* Process T2 and get Q2 and λ2. */
    matrix_divide_and_conquer_inner(mat, m, e_ind, &Q2, &lambda2, eps);

    /* Calculate P, P_t and D. */
    D = matrix_from_two_blocks(lambda1, lambda2);
    P = matrix_get_permut(D);    /* P.    */
    P_t = matrix_transpose(P);   /* P_t.  */
    tmp = matrix_mul(P, D);
    matrix_free(D);
    D = tmp;
    tmp = matrix_mul(D, P_t);
    matrix_free(D);
    D = tmp;                     /* D.    */

    /* Calculate Q, Q_t and u. */
    (*Q) = matrix_from_two_blocks(Q1, Q2);
    Q_t = matrix_transpose(*Q);
    {
        v = matrix_create(matrix_height(mat), 1, 0.0);
        matrix_set(v, m - 1, 0, 1.0);
        matrix_set(v, m,     0, 1.0);
    }
    u = matrix_mul(Q_t, v);

    tmp = matrix_mul((*Q), P_t);
    matrix_free(*Q);
    (*Q) = tmp;                  /* Q.    */

    matrix_free(Q_t);
    Q_t = matrix_transpose(*Q);  /* Q_t.  */

    tmp = matrix_mul(P, u);
    matrix_free(u);
    u = tmp;                     /* u.    */

    /* Calculate λ. */
    lambdas = (double*) malloc(matrix_height(D) * sizeof(double));
    printf("-----------------------------\n");
    printf("D:\n");
    matrix_print(D);
    printf("u:\n");
    matrix_print(u);
    printf("lambdas:\n");
    for (i = 0; i < matrix_height(D) - 1; i++) {
        lambdas[i] = solve_secular_equation_common(D, u, b_m, i, eps);
        printf("%lf ", lambdas[i]);
    }
    if (b_m > 0.0) {
        lambdas[matrix_height(D) - 1] = solve_secular_equation_d0_plus_inf(D, u, b_m, eps);
        printf("(b_m > 0) %lf ", lambdas[matrix_height(D) - 1]);
    } else {
        lambdas[matrix_height(D) - 1] = solve_secular_equation_minus_inf_dn(D, u, b_m, eps);
        printf("(b_m < 0) %lf ", lambdas[matrix_height(D) - 1]);
    }
    printf("\n");

    printf("-----------------------------\n");

    /* Restore original matrix from T1. */
    matrix_set(mat, m - 1, m - 1, matrix_get(mat, m - 1, m - 1) + b_m);
    /* Restore original matrix from T2. */
    matrix_set(mat, m, m, matrix_get(mat, m, m) + b_m);

    free(lambdas);

    matrix_free(Q1);
    matrix_free(lambda1);

    matrix_free(Q2);
    matrix_free(lambda2);
}

void matrix_divide_and_conquer(matrix_type_t mat, double eps)
{
    matrix_type_t Q;
    matrix_type_t lambda;

    matrix_divide_and_conquer_inner(mat, 0, matrix_height(mat), &Q, &lambda, eps);

    printf("Q:\n");
    matrix_print(Q);
    printf("lambda:\n");
    matrix_print(lambda);
}
