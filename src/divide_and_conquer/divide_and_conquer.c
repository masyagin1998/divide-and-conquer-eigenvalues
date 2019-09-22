#include "divide_and_conquer.h"

#include "matrix_algo.h"
#include "secular_equation.h"
#include "macro.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define DEBUG_STEPS

#ifdef DEBUG_STEPS
#include "utils.h"
#endif  /* DEBUG_STEPS */

enum LAMBDA_TYPE
{
    LAMBDA_TYPE_NORMAL,          /* lambda, found by solving secular equation for (d_i, d_i+1). */
    LAMBDA_TYPE_DEFLATED_ZERO_U, /* deflated lambda with u[i][0] = 0.                           */
    LAMBDA_TYPE_DEFLATED_EQ_D,   /* deflated lambda with d[i][i] = d[i+1][i+1].                 */
};

struct LAMBDA
{
    int              idx;
    enum LAMBDA_TYPE type;
    long double           lambda;
};

struct D_I_J
{
    unsigned j;
    long double   d_i_j;
    unsigned k_d_i_j;
    unsigned l_j;
};

int exists_in_d_i_j(const struct D_I_J*d_i_j_s, unsigned d_i_j_s_len, long double lambda, long double eps)
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

int is_d_i_equal_to_d_i_j(const struct D_I_J*d_i_j_s, unsigned d_i_j_s_len, long double d_i, long double eps)
{
    unsigned k;

    for (k = 0; k < d_i_j_s_len; k++) {
        if (are_equal(d_i_j_s[k].d_i_j, d_i)) {
            return k;
        }
    }

    return -1;
}

static long double calculate_w_i_0(const matrix_type_t D, const struct LAMBDA*lambdas, unsigned i)
{
    unsigned j;
    
    long double num = 1.0;
    long double denom = 1.0;
    long double res;

    for (j = 0; j < matrix_height(D); j++) {
        printf("lambda: %Lf; d: %Lf\n", lambdas[j].lambda, matrix_get(D, i, i));
        num *= (lambdas[j].lambda - matrix_get(D, i, i));
        if (j != i) {
            printf("d_j: %Lf, d_i: %Lf\n", matrix_get(D, j, j), matrix_get(D, i, i));
            denom *= (matrix_get(D, j, j) - matrix_get(D, i, i));
        }
    }

    printf("num: %Lf; denom: %Lf\n", num, denom);

    res = num / denom;
    if (res < 0.0) {
        res = -res;
    }
    res = sqrt(res);

    printf("res: %Lf\n", res);

    return res;
}

static long double calculate_w_i_1(const matrix_type_t D, const struct LAMBDA*lambdas, unsigned i, long double k_d_i_j, long double eps)
{
    unsigned j;

    long double num = 1.0;
    long double denom = 1.0;
    long double res;

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
    if (res < 0) {
        res = -res;
    }    
    res = sqrt(res);

    return res;
    
}

static matrix_type_t matrix_from_eigs(const matrix_type_t*eigenvectors, unsigned len)
{
    unsigned i, j;
    
    matrix_type_t res = matrix_create(len, len, 0.0);
    if (res == NULL) {
        goto err0;
    }

    for (i = 0; i < len; i++) {
        for (j = 0; j < len; j++) {
            long double cell = matrix_get(eigenvectors[i], j, 0);
            matrix_set(res, j, i, cell);
        }
    }

    return res;

 err0:
    return NULL;
}

static matrix_type_t matrix_from_lambdas(const struct LAMBDA*lambdas, unsigned len)
{
    unsigned i;

    matrix_type_t res = matrix_create(len, len, 0.0);
    if (res == NULL) {
        goto err0;
    }

    for (i = 0; i < len; i++) {
        long double cell = lambdas[i].lambda;
        matrix_set(res, i, i, cell);
    }

    return res;

 err0:
    return NULL;    
}

static void matrix_divide_and_conquer_inner(matrix_type_t mat, matrix_type_t*Q, matrix_type_t*LAMBDA, long double eps, unsigned step)
{    
    unsigned i, j;
    
    unsigned m;
    long double b_m;

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
    
    if (matrix_height(mat) == 1) {
        (*LAMBDA) = matrix_create(1, 1, matrix_get(mat, 0, 0));
        (*Q) = matrix_create(1, 1, 1.0);
        return;
    }

    /* m and b_m. */
    m = matrix_height(mat) / 2;
    b_m = matrix_get(mat, m, m - 1);

    /* Create T1, process it, and get Q1 and λ1. */
    matrix_type_t mat1 = matrix_create(m, m, 0.0);
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            long double cell = matrix_get(mat, i, j);
            matrix_set(mat1, i, j, cell);
        }
    }
    matrix_set(mat1, m - 1, m - 1, matrix_get(mat1, m - 1, m - 1) - b_m);    
    matrix_divide_and_conquer_inner(mat1, &Q1, &lambda1, eps, step * 2);
    /* Create T2, process it and get Q2 and λ2. */
    matrix_type_t mat2 = matrix_create(matrix_height(mat) - m, matrix_width(mat) - m, 0.0);
    for (i = m; i < matrix_height(mat); i++) {
        for (j = m; j < matrix_width(mat); j++) {
            long double cell = matrix_get(mat, i, j);
            matrix_set(mat2, i - m, j - m, cell);
        }
    }
    matrix_set(mat2, 0, 0, matrix_get(mat2, 0, 0) - b_m);
    matrix_divide_and_conquer_inner(mat2, &Q2, &lambda2, eps, step * 2 + 1);

#ifdef DEBUG_STEPS
    printf("Step: %u\n", step);
    printf("In:\n");
    matrix_print(mat);
    printf("B_m:\n");
    printf("%Lf\n", b_m);    
#endif  /* DEBUG_STEPS */

    /* Calculate P, P_t and D. */
    D = matrix_from_two_blocks(lambda1, lambda2);
#ifdef DEBUG_STEPS
    printf("D:\n");
    matrix_print(D);
#endif  /* DEBUG_STEPS */
    P = matrix_diag_permut(D);   /* P.    */
    P_t = matrix_transpose(P);   /* P_t.  */
    tmp = matrix_mul(P, D);
    matrix_free(D);
    D = tmp;
    tmp = matrix_mul(D, P_t);
    matrix_free(D);
    D = tmp;                     /* D.    */

#ifdef DEBUG_STEPS
    printf("P:\n");
    matrix_print(P);
    printf("P_t:\n");
    matrix_print(P_t);
    printf("D_new:\n");
    matrix_print(D);
#endif  /* DEBUG_STEPS */

    /* Calculate Q, Q_t and u. */
    (*Q) = matrix_from_two_blocks(Q1, Q2);
#ifdef DEBUG_STEPS    
    printf("Q:\n");
    matrix_print(*Q);
#endif  /* DEBUG_STEPS */
    Q_t = matrix_transpose(*Q);
#ifdef DEBUG_STEPS
    printf("Q_t:\n");
    matrix_print(Q_t);
#endif  /* DEBUG_STEPS */
    
    v = matrix_create(matrix_height(mat), 1, 0.0);
    matrix_set(v, m - 1, 0, 1.0);
    matrix_set(v, m,     0, 1.0);

#ifdef DEBUG_STEPS
    printf("V:\n");
    matrix_print(v);
#endif  /* DEBUG_STEPS */
    
    u = matrix_mul(Q_t, v);
#ifdef DEBUG_STEPS    
    printf("U:\n");
    matrix_print(u);
#endif  /* DEBUG_STEPS */    
    matrix_free(Q_t);
    matrix_free(v);

    tmp = matrix_mul((*Q), P_t);
    matrix_free(*Q);
    matrix_free(P_t);
    (*Q) = tmp;                  /* Q.    */
#ifdef DEBUG_STEPS    
    printf("Q_new:\n");
    matrix_print(*Q);
#endif  /* DEBUG_STEPS */    

    tmp = matrix_mul(P, u);
    matrix_free(P);
    matrix_free(u);
    u = tmp;                     /* u.    */
#ifdef DEBUG_STEPS
    printf("U new:\n");
    matrix_print(u);
#endif  /* DEBUG_STEPS */    

    /* Calculate λ. */
#ifdef DEBUG_STEPS    
    printf("Lambdas:\n");
#endif  /* DEBUG_STEPS */
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
#ifdef DEBUG_STEPS
        printf("(d0, +oo): d0 = %Lf; lambda = %Lf; type = %u;\n",
               matrix_get(D, 0, 0),
               lambdas[lambdas_len].lambda, lambdas[lambdas_len].type);
#endif  /* DEBUG_STEPS */
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
#ifdef DEBUG_STEPS
        printf("(d%u, d%u): d%u = %Lf; d%u = %Lf; lambda = %Lf; type = %u;\n",
               i + 1, i, i + 1, matrix_get(D, i + 1, i + 1), i, matrix_get(D, i, i),
               lambdas[lambdas_len].lambda, lambdas[lambdas_len].type);
#endif  /* DEBUG_STEPS */
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
#ifdef DEBUG_STEPS
        printf("(-oo, d_%u): d_%u = %Lf; lambda = %Lf; type = %u;\n",
               matrix_height(D) - 1, matrix_height(D) - 1,                   
               matrix_get(D, matrix_height(D) - 1, matrix_height(D) - 1),
               lambdas[lambdas_len].lambda, lambdas[lambdas_len].type);
#endif  /* DEBUG_STEPS */
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
            long double cell = calculate_w_i_0(D, lambdas, i);
            matrix_set(w, i, 0, cell);
        } else if ((idx_d_i_j != -1) && (d_i_j_s[idx_d_i_j].l_j == d_i_j_s[idx_d_i_j].k_d_i_j - 1) && (d_i_j_s[idx_d_i_j].k_d_i_j - 1 > 0)) {
            long double cell = calculate_w_i_1(D, lambdas, i, d_i_j_s[idx_d_i_j].k_d_i_j, eps);
            matrix_set(w, i, 0, cell);
        } else {
            printf("\nROR\n");
            exit(1);
        }
    }

#ifdef DEBUG_STEPS
    {
        matrix_type_t w_t;
        matrix_type_t w_x_w_t;
        matrix_type_t D_circum;
        
        printf("W:\n");
        matrix_print(w);
        w_t = matrix_transpose(w);
        printf("W_t:\n");
        matrix_print(w_t);
        w_x_w_t = matrix_mul(w, w_t);
        matrix_free(w_t);
        printf("W x W_t:\n");
        matrix_print(w_x_w_t);
        if (b_m > 0.0f) {
            D_circum = matrix_plus(D, w_x_w_t);
        } else {
            D_circum = matrix_minus(D, w_x_w_t);
        }
        matrix_free(w_x_w_t);
        printf("D_circum:\n");
        matrix_print(D_circum);
        matrix_free(D_circum);
    }
#endif  /* DEBUG_STEPS */

    matrix_type_t evs[matrix_height(D)];
    for (i = 0; i < lambdas_len; i++) {
        matrix_type_t a_i_E;
        matrix_type_t m;
        matrix_type_t m_inv;
        
        a_i_E = matrix_diag(matrix_height(D), lambdas[i].lambda);
        m     = matrix_minus(a_i_E, D);
#ifdef DEBUG_STEPS
        printf("D - a_i * E:\n");
        matrix_print(m);
#endif  /* DEBUG_STEPS */        
        matrix_free(a_i_E);
        m_inv = matrix_inverse(m);
#ifdef DEBUG_STEPS
        printf("(D - a_i * E)^(-1):\n");
        matrix_print(m_inv);
#endif  /* DEBUG_STEPS */        
        matrix_free(m);
        evs[i] = matrix_mul(m_inv, w);
        long double len = vector_len(evs[i]);
        for (j = 0; j < matrix_height(evs[i]); j++) {
            long double cell = matrix_get(evs[i], j, 0) / len;
            matrix_set(evs[i], j, 0, cell);
        }
        matrix_free(m_inv);
#ifdef DEBUG_STEPS
        printf("eigenvector:\n");
        matrix_print(evs[i]);
#endif  /* DEBUG_STEPS */
    }
    
    matrix_free(w);
    
    matrix_type_t eigs = matrix_from_eigs(evs, matrix_height(D));

    for (i = 0; i < matrix_height(D); i++) {
        matrix_free(evs[i]);
    }

    (*LAMBDA) = matrix_from_lambdas(lambdas, lambdas_len);

#ifdef DEBUG_STEPS
        printf("eigenvectors:\n");
        matrix_print(eigs);
        printf("lambdas:\n");
        matrix_print(*LAMBDA);
        printf("eigs * lambdas * eigs_t\n");
        matrix_print(matrix_mul(matrix_mul(eigs, *LAMBDA), matrix_transpose(eigs)));
#endif  /* DEBUG_STEPS */    
    
    tmp = matrix_mul((*Q), eigs);
    matrix_free(*Q);
    matrix_free(eigs);
    (*Q) = tmp;

    matrix_free(D);

    free(lambdas);
    free(d_i_j_s);

    matrix_free(Q1);
    matrix_free(lambda1);

    matrix_free(Q2);
    matrix_free(lambda2);
}

void matrix_divide_and_conquer(matrix_type_t mat, long double eps)
{
    matrix_type_t Q;
    matrix_type_t lambda;
    matrix_type_t Q_t;
    matrix_type_t res;
    matrix_type_t tmp;

    matrix_divide_and_conquer_inner(mat, &Q, &lambda, eps, 1);
    Q_t = matrix_inverse(Q);
    tmp = matrix_mul(Q, lambda);
    res = matrix_mul(tmp, Q_t);
    
    matrix_free(tmp);

#ifdef DEBUG_STEPS    
    printf("\n");
    printf("Q:\n");
    matrix_print(Q);
    printf("LAMBDA:\n");
    matrix_print(lambda);
    printf("Q_t:\n");
    matrix_print(Q_t);
    printf("RES:\n");
    matrix_print(res);
#endif  /* DEBUG_STEPS */

    matrix_free(Q);
    matrix_free(lambda);
    matrix_free(Q_t);
    matrix_free(res);
}
