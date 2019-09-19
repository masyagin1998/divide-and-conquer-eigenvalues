#include "divide_and_conquer.h"
#include "secular_equation.h"
#include "macro.h"
#include "utils.h"

#include <math.h>
#include <stdlib.h>

static matrix_type_t matrix_from_two_blocks(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(matrix_height(a) + matrix_height(b), matrix_width(a) + matrix_width(b), 0.0);
    if (res == NULL) {
        goto err0;
    }
    
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

 err0:
    return NULL;
}

static matrix_type_t matrix_transpose(const matrix_type_t mat)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(matrix_width(mat), matrix_height(mat), 0.0);
    if (res == NULL) {
        goto err0;
    }

    for (i = 0; i < matrix_width(mat); i++) {
        for (j = 0; j < matrix_height(mat); j++) {
            matrix_set(res, i, j, matrix_get(mat, j, i));
        }
    }

    return res;

 err0:
    return NULL;
}

static matrix_type_t matrix_mul(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j, k;
    
    matrix_type_t res = matrix_create(matrix_height(a), matrix_width(b), 0.0);
    if (res == NULL) {
        goto err0;
    }
    
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

 err0:
    return NULL;
}

static matrix_type_t matrix_permut(const matrix_type_t mat)
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

static matrix_type_t matrix_plus(const matrix_type_t a, const matrix_type_t b)
{
    unsigned i, j;

    matrix_type_t res = matrix_create(matrix_height(a), matrix_width(a), 0.0);
    if (res == NULL) {
        goto err0;
    }

    for (i = 0; i < matrix_height(res); i++) {
        for (j = 0; j < matrix_width(res); j++) {
            double cell = matrix_get(a, i, j) + matrix_get(b, i, j);
            matrix_set(res, i, j, cell);
        }
    }

    return res;

 err0:
    return NULL;
}

enum LAMBDA_TYPE
{
    LAMBDA_TYPE_NORMAL,          /* lambda, found by solving secular equation for (d_i, d_i+1). */
    LAMBDA_TYPE_DEFLATED_ZERO_U, /* deflated lambda with u[i][0] = 0.                           */
    LAMBDA_TYPE_DEFLATED_EQ_D,   /* deflated lambda with d[i][i] = d[i+1][i+1].                 */
};

#include <stdio.h>

struct LAMBDA
{
    int              idx;
    enum LAMBDA_TYPE type;
    double           lambda;
};

struct D_I_J
{
    unsigned j;
    double   d_i_j;
    unsigned k_d_i_j;
    unsigned l_j;
};

int exists_in_d_i_j(const struct D_I_J*d_i_j_s, unsigned d_i_j_s_len, double lambda, double eps)
{
    unsigned i;
    
    for (i = 0; i < d_i_j_s_len; i++) {
        if (are_equal(d_i_j_s[i].d_i_j, lambda)) {
            return i;
        }
    }

    return -1;
}

int is_i_equal_to_i_j(const struct D_I_J*d_i_j_s, unsigned d_i_j_s_len, unsigned i)
{
    unsigned k;
    
    for (k = 0; k < d_i_j_s_len; k++) {
        if (d_i_j_s[k].j == i) {
            return k;
        }
    }

    return -1;
}

int is_d_i_equal_to_d_i_j(const struct D_I_J*d_i_j_s, unsigned d_i_j_s_len, double d_i, double eps)
{
    unsigned k;

    for (k = 0; k < d_i_j_s_len; k++) {
        if (are_equal(d_i_j_s[k].d_i_j, d_i)) {
            return k;
        }
    }

    return -1;
}

static double calculate_w_i_0(const matrix_type_t D, const struct LAMBDA*lambdas, unsigned i)
{
    unsigned j;
    
    double num = 1.0;
    double denom = 1.0;
    double res;

    for (j = 0; j < matrix_height(D); j++) {
        num *= (lambdas[j].lambda - matrix_get(D, i, i));
        if (j != i) {
            denom *= (matrix_get(D, j, j) - matrix_get(D, i, i));
        }
    }

    res = num / denom;
    res = sqrt(res);

    return res;
}

static double calculate_w_i_1(const matrix_type_t D, const struct LAMBDA*lambdas, unsigned i, double k_d_i_j, double eps)
{
    unsigned j;

    double num = 1.0;
    double denom = 1.0;
    double res;

    for (j = 0; j < matrix_height(D); j++) {
        if (!are_equal(lambdas[j].lambda, matrix_get(D, i, i))) {
            num *= (lambdas[j].lambda - matrix_get(D, i, i));
        }
        if (!are_equal(matrix_get(D, j, j), matrix_get(D, i, i))) {
            denom *= (matrix_get(D, j, j) - matrix_get(D, i, i));
        }
    }

    denom *= k_d_i_j;

    res = num / denom;
    res = sqrt(res);

    return res;
    
}

static void matrix_divide_and_conquer_inner(matrix_type_t mat, unsigned b_ind, unsigned e_ind, matrix_type_t*Q, matrix_type_t*LAMBDA, double eps)
{
    unsigned i, j;
    
    unsigned m;
    double b_m;

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

    struct LAMBDA*lambdas;
    unsigned      lambdas_len;

    struct D_I_J*d_i_j_s;
    unsigned     d_i_j_s_len;

    matrix_type_t w;
    matrix_type_t w_t;
    matrix_type_t w_x_w_t;
    
    if ((e_ind - b_ind == 1)) {
        (*LAMBDA) = matrix_create(1, 1, matrix_get(mat, b_ind, b_ind));
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
    P = matrix_permut(D);        /* P.    */
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
    
    v = matrix_create(matrix_height(mat), 1, 0.0);
    matrix_set(v, m - 1, 0, 1.0);
    matrix_set(v, m,     0, 1.0);
    u = matrix_mul(Q_t, v);
    matrix_free(Q_t);
    matrix_free(v);

    tmp = matrix_mul((*Q), P_t);
    matrix_free(*Q);
    matrix_free(P_t);
    (*Q) = tmp;                  /* Q.    */

    tmp = matrix_mul(P, u);
    matrix_free(P);
    matrix_free(u);
    u = tmp;                     /* u.    */

    /* Calculate λ. */
    lambdas = (struct LAMBDA*) malloc(matrix_height(D) * sizeof(struct LAMBDA));
    lambdas_len = 0;

    if (b_m > 0.0) {
        lambdas[lambdas_len].idx = -1;
        if (is_zero(matrix_get(u, 0, 0))) {
            lambdas[lambdas_len].type = LAMBDA_TYPE_DEFLATED_ZERO_U;
            lambdas[lambdas_len].lambda = matrix_get(D, 0, 0);
        } else {
            lambdas[lambdas_len].type = LAMBDA_TYPE_NORMAL;
            lambdas[lambdas_len].lambda = solve_secular_equation_d0_plus_inf(D, u, b_m, eps);
        }
        lambdas_len++;
    }

    for (i = 0; i < matrix_height(D) - 1; i++) {
        lambdas[lambdas_len].idx = i;
        if (is_zero(matrix_get(u, i, 0))) {
            lambdas[lambdas_len].type = LAMBDA_TYPE_DEFLATED_ZERO_U;
            lambdas[lambdas_len].lambda = matrix_get(D, i, i);
        } else if (are_equal(matrix_get(D, i, i), matrix_get(D, i + 1, i + 1))) {
            lambdas[lambdas_len].type = LAMBDA_TYPE_DEFLATED_EQ_D;
            lambdas[lambdas_len].lambda = matrix_get(D, i, i);
        } else {
            lambdas[lambdas_len].type = LAMBDA_TYPE_NORMAL;
            lambdas[lambdas_len].lambda = solve_secular_equation_common(D, u, b_m, i, eps);
        }
        lambdas_len++;
    }

    if (b_m < 0.0) {
        lambdas[lambdas_len].idx = matrix_height(D) - 1;
        if (is_zero(matrix_get(u, matrix_height(D) - 1, 0))) {
            lambdas[lambdas_len].type = LAMBDA_TYPE_DEFLATED_ZERO_U;
            lambdas[lambdas_len].lambda = matrix_get(D, matrix_height(D) - 1, matrix_height(D) - 1);
        } else {
            lambdas[lambdas_len].type = LAMBDA_TYPE_NORMAL;
            lambdas[lambdas_len].lambda = solve_secular_equation_minus_inf_dn(D, u, b_m, eps);
        }
        lambdas_len++;
    }

    matrix_free(u);

    /* Calculate D_i_j. */
    d_i_j_s = (struct D_I_J*) malloc(matrix_height(D) * sizeof(struct D_I_J));
    d_i_j_s_len = 0;

    for (i = 0; i < matrix_height(D); i++) {
        int idx = exists_in_d_i_j(d_i_j_s, d_i_j_s_len, matrix_get(D, i, i), eps);
        if (idx == -1) {
            d_i_j_s[d_i_j_s_len].j = i;
            d_i_j_s[d_i_j_s_len].d_i_j = matrix_get(D, i, i);
            d_i_j_s[d_i_j_s_len].k_d_i_j = 1;
            d_i_j_s[d_i_j_s_len].l_j = 0;
            for (j = 0; j < lambdas_len; j++) {
                if (are_equal(lambdas[j].lambda, d_i_j_s[d_i_j_s_len].d_i_j)) {
                    d_i_j_s[d_i_j_s_len].l_j++;
                }
            }
            d_i_j_s_len++;
        } else {
            d_i_j_s[idx].k_d_i_j++;
        }
    }

    /* Calculate w. */
    w = matrix_create(matrix_height(D), 1, 0.0);

    for (i = 0; i < matrix_height(w); i++) {
        int idx_i_j = is_i_equal_to_i_j(d_i_j_s, d_i_j_s_len, i);
        int idx_d_i_j = is_d_i_equal_to_d_i_j(d_i_j_s, d_i_j_s_len, matrix_get(D, i, i), eps);
        if (((idx_i_j != -1) && (d_i_j_s[idx_i_j].l_j >= d_i_j_s[idx_i_j].k_d_i_j) && (d_i_j_s[idx_i_j].k_d_i_j == 1)) ||
            ((idx_d_i_j != -1) && (d_i_j_s[idx_d_i_j].l_j > d_i_j_s[idx_d_i_j].k_d_i_j - 1) && (d_i_j_s[idx_d_i_j].k_d_i_j - 1 > 0))) {
            matrix_set(w, i, 0, 0.0);
        } else if ((idx_i_j != -1) && (d_i_j_s[idx_i_j].l_j == 0) && (d_i_j_s[idx_i_j].k_d_i_j == 1)) {
            double cell = calculate_w_i_0(D, lambdas, i);
            matrix_set(w, i, 0, cell);
        } else if ((idx_d_i_j != -1) && (d_i_j_s[idx_d_i_j].l_j == d_i_j_s[idx_d_i_j].k_d_i_j - 1) && (d_i_j_s[idx_d_i_j].k_d_i_j - 1 > 0)) {
            double cell = calculate_w_i_1(D, lambdas, i, d_i_j_s[idx_d_i_j].k_d_i_j, eps);
            matrix_set(w, i, 0, cell);            
        } else {
            printf("\nROR\n");
            exit(1);
        }
    }

    w_t = matrix_transpose(w);
    w_x_w_t = matrix_mul(w, w_t);
    matrix_free(w);
    matrix_free(w_t);
    (*LAMBDA) = matrix_plus(D, w_x_w_t);
    matrix_free(w_x_w_t);
    matrix_free(D);

    /* Restore original matrix from T1. */
    matrix_set(mat, m - 1, m - 1, matrix_get(mat, m - 1, m - 1) + b_m);
    /* Restore original matrix from T2. */
    matrix_set(mat, m, m, matrix_get(mat, m, m) + b_m);

    free(lambdas);
    free(d_i_j_s);

    matrix_free(Q1);
    matrix_free(lambda1);

    matrix_free(Q2);
    matrix_free(lambda2);
}

void matrix_divide_and_conquer(matrix_type_t mat, double eps)
{
    matrix_type_t Q;
    matrix_type_t lambda;
    matrix_type_t Q_t;
    matrix_type_t res;
    matrix_type_t tmp;

    matrix_divide_and_conquer_inner(mat, 0, matrix_height(mat), &Q, &lambda, eps);
    Q_t = matrix_transpose(Q);
    tmp = matrix_mul(Q, lambda);
    res = matrix_mul(tmp, Q_t);
    matrix_free(tmp);

    printf("\n");
    printf("Q:\n");
    matrix_print(Q);
    printf("LAMBDA:\n");
    matrix_print(lambda);
    printf("Q_T:\n");
    matrix_print(Q_t);
    printf("RES:\n");
    matrix_print(res);

    matrix_free(Q);
    matrix_free(lambda);
    matrix_free(Q_t);
    matrix_free(res);
}
