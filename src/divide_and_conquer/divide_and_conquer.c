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

enum LAMBDA_TYPE
{
    LAMBDA_TYPE_NORMAL,          /* lambda, found by solving secular equation for (d_i, d_i+1). */
    LAMBDA_TYPE_DEFLATED_ZERO_U, /* deflated lambda with u[i][0] = 0.                           */
    LAMBDA_TYPE_DEFLATED_EQ_D,   /* deflated lambda with d[i][i] = d[i+1][i+1].                 */
};

#define is_zero(v) ((v) < eps)
#define are_equal(v1, v2) ((fabs((v1) - (v2))) < eps)

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

    enum LAMBDA_TYPE *lambda_types;
    double           *lambdas;
    unsigned         lambdas_len;
    
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
    lambda_types = (enum LAMBDA_TYPE*) malloc(matrix_height(D) * sizeof(enum LAMBDA_TYPE));
    lambdas = (double*) malloc(matrix_height(D) * sizeof(double));
    lambdas_len = 0;
    printf("-----------------------------\n");
    printf("b_m:\n");
    printf("%lf\n", b_m);
    printf("D:\n");
    matrix_print(D);
    printf("u:\n");
    matrix_print(u);

    for (i = 0; i < matrix_height(D) - 1; i++) {
        if (is_zero(matrix_get(u, i, 0))) {
            lambda_types[lambdas_len] = LAMBDA_TYPE_DEFLATED_ZERO_U;
            lambdas[lambdas_len] = matrix_get(D, i, i);
        } else if (are_equal(matrix_get(D, i, i), matrix_get(D, i + 1, i + 1))) {
            lambda_types[lambdas_len] = LAMBDA_TYPE_DEFLATED_EQ_D;
            lambdas[lambdas_len] = matrix_get(D, i, i);
        } else {
            lambda_types[lambdas_len] = LAMBDA_TYPE_NORMAL;
            lambdas[lambdas_len] = solve_secular_equation_common(D, u, b_m, i, eps);
        }
        lambdas_len++;
    }
    if (is_zero(matrix_get(u, matrix_height(D) - 1, 0))) {
        lambda_types[lambdas_len] = LAMBDA_TYPE_DEFLATED_ZERO_U;
        lambdas[lambdas_len] = matrix_get(D, matrix_height(D) - 1, matrix_height(D) - 1);
    } else {
        if (b_m > 0.0) {
            lambda_types[lambdas_len] = LAMBDA_TYPE_NORMAL;
            lambdas[lambdas_len] = solve_secular_equation_d0_plus_inf(D, u, b_m, eps);
        } else {
            lambda_types[lambdas_len] = LAMBDA_TYPE_NORMAL;
            lambdas[lambdas_len] = solve_secular_equation_minus_inf_dn(D, u, b_m, eps);
        }
    }
    lambdas_len++;

    printf("lambdas:\n");
    for (i = 0; i < lambdas_len; i++) {
        printf("lambda[%d] = %lf (type: %u); ", i, lambdas[i], lambda_types[i]);
    }
    printf("\n");

    printf("-----------------------------\n");

    /* */

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
