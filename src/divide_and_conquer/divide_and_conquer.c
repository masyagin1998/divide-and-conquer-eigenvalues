#include "divide_and_conquer.h"

#include "macro.h"
#include "matrix_algo.h"
#include "deflation.h"
#include "secular_equation.h"

#include <stddef.h>
#include <assert.h>
#include <math.h>

/*
  If DIVIDE_AND_CONQUER_DEBUG is defined,
  matrix_divide_and_conquer will print
  some debug info to stdout.
 */
/* #define DIVIDE_AND_CONQUER_DEBUG */

#ifdef DIVIDE_AND_CONQUER_DEBUG
#include <stdio.h>
#include "utils.h"
#endif  /* DIVIDE_AND_CONQUER_DEBUG */

static matrix_type_t matrix_from_part(const matrix_type_t T, unsigned b, unsigned e)
{
    unsigned n = (e - b + 1);
    unsigned i, j;

    matrix_type_t res = matrix_create(n, n);
    if (res == NULL) {
        goto err0;
    }
    
    for (i = b; i <= e; i++) {
        for (j = b; j <= e; j++) {
            long double cell = matrix_get(T, i, j);
            matrix_set(res, i - b + 1, j - b + 1, cell);
        }
    }

    return res;

 err0:
    return NULL;
}

int matrix_divide_and_conquer(const matrix_type_t T, matrix_type_t*Q, matrix_type_t*L, long double eps)
{
    int r;
    
    unsigned n, m;
    unsigned i, j;
    
    long double rho;
    long double cell;
    matrix_type_t T1, Q1, L1, T2, Q2, L2;
    matrix_type_t D;
    int sign;
    matrix_type_t v;
    matrix_type_t P, P1;
    matrix_type_t tmp1, tmp2;
    matrix_type_t v_prime;
    matrix_type_t eigenvalues, eigenvectors;
    unsigned n_deflated;
    matrix_type_t G;
    matrix_type_t lambda;
    matrix_type_t v_hat;
    unsigned k;
    unsigned jj;
    
    assert(matrix_height(T) == matrix_width(T));

    n = matrix_height(T);
    
    if (n == 1) {
        (*Q) = matrix_create(1, 1);
        if ((*Q) == NULL) {
            return -1;
        }
        matrix_set((*Q), 1, 1, 1.0);
        (*L) = matrix_create(1, 1);
        if ((*L) == NULL) {
            matrix_free((*Q));
            return -2;
        }        
        matrix_set((*L), 1, 1, matrix_get(T, 1, 1));
        return 0;
    }

    if (n % 2 == 0) {
        m = n / 2;
    } else {
        m = (n - 1) / 2;
    }

    rho = fabsl(matrix_get(T, m, m + 1));

    T1 = matrix_from_part(T, 1, m);
    if (T1 == NULL) {
        r = -3;
        goto err0;
    }
    cell = matrix_get(T1, m, m);
    matrix_set(T1, m, m, cell - rho);

    T2 = matrix_from_part(T, m + 1, n);
    if (T2 == NULL) {
        r = -4;
        goto err1;
    }
    cell = matrix_get(T2, 1, 1);
    matrix_set(T2, 1, 1, cell - rho);

    r = matrix_divide_and_conquer(T1, &Q1, &L1, eps);
    if (r) {
        r = -5;
        goto err1;
    }
    r = matrix_divide_and_conquer(T2, &Q2, &L2, eps);
    if (r) {
        r = -6;
        goto err1;
    }

    D = matrix_from_two_blocks(L1, L2);
    if (D == NULL) {
        r = -7;
        goto err1;
    }

    v = matrix_create(n, 1);
    if (v == NULL) {
        r = -8;
        goto err2;
    }
    sign = (matrix_get(T, m, m + 1)) > 0.0 ? 1 : -1;
    for (i = 1; i <= m; i++) {
        cell = matrix_get(Q1, m, i) * sign;
        matrix_set(v, i, 1, cell);
    }
    for (i = m + 1; i <= n; i++) {
        cell = matrix_get(Q2, 1, i - m);
        matrix_set(v, i, 1, cell);
    }

#ifdef DIVIDE_AND_CONQUER_DEBUG
    printf("----divide-and-conquer-beg----\n");
    printf("T1:\n");
    matrix_print(T1);
    printf("Q1:\n");
    matrix_print(Q1);
    printf("L1:\n");
    matrix_print(L1);
    printf("\n");
    printf("T2:\n");
    matrix_print(T2);
    printf("Q2:\n");
    matrix_print(Q2);
    printf("L2:\n");
    matrix_print(L2);
    printf("\n");
    printf("D:\n");
    matrix_print(D);
    printf("v:\n");
    matrix_print(v);
    printf("\n");
#endif  /* DIVIDE_AND_CONQUER_DEBUG */

    P = matrix_diag_permut(D);
    tmp1 = matrix_mul(P, D);
    tmp2 = matrix_transpose(P);
    matrix_free(D);
    D = matrix_mul(tmp1, tmp2);
    matrix_free(tmp1);
    matrix_free(tmp2);

    tmp1 = matrix_mul(P, v);
    matrix_free(v);
    v = tmp1;

#ifdef DIVIDE_AND_CONQUER_DEBUG
    printf("P:\n");
    matrix_print(P);
    printf("D new:\n");
    matrix_print(D);
    printf("v new:\n");
    matrix_print(v);
    printf("\n");
#endif  /* DIVIDE_AND_CONQUER_DEBUG */

    r = deflate(&D, &v, &v_prime, &eigenvalues, &eigenvectors, &n_deflated, &G, eps);
    if (r) {
        r = -9;
        goto err2;
    }

#ifdef DIVIDE_AND_CONQUER_DEBUG
    printf("v prime deflate:\n");
    matrix_print(v_prime);
    printf("eigenvalues deflate:\n");
    matrix_print(eigenvalues);
    printf("eigenvectors deflate:\n");
    matrix_print(eigenvectors);
    printf("n deflated:\n");
    printf("%u\n", n_deflated);
    printf("G:\n");
    matrix_print(G);
#endif  /* DIVIDE_AND_CONQUER_DEBUG */

    lambda = matrix_create(n_deflated, 1);
    if (lambda == NULL) {
        r = -10;
        goto err3;
    }

    tmp1 = matrix_transpose(v_prime);

    for (i = 1; i <= n_deflated; i++) {
        long double dn;
        long double lambda_init;
        long double cell;
        if (i == n_deflated) {
            long double vn = vector_norm(v);
            dn = matrix_get(D, n_deflated, n_deflated) + vn * vn;
            lambda_init = (matrix_get(D, n_deflated, n_deflated) + dn) / 2.0;
        } else {
            lambda_init = (matrix_get(D, i, i) + matrix_get(D, i + 1, i + 1)) / 2.0;
        }
        cell = solve_secular_equation(rho, D, tmp1, i, lambda_init, n_deflated, eps);
        matrix_set(lambda, i, 1, cell);
    }

    matrix_free(tmp1);

#ifdef DIVIDE_AND_CONQUER_DEBUG
    printf("v_prime:\n");
    matrix_print(v_prime);
    printf("lambda:\n");
    matrix_print(lambda);
#endif  /* DIVIDE_AND_CONQUER_DEBUG */

    v_hat = matrix_create(n_deflated, 1);
    for (k = 1; k <= n_deflated; k++) {
        long double cell;
        long double a;
        long double b;
        int sign = matrix_get(v_prime, 1, k) > 0.0 ? 1 : -1;
        if (k == 1) {
            a = 1.0;
            for (i = 1; i <= matrix_height(lambda); i++) {
                a *= (matrix_get(lambda, i, 1) - matrix_get(D, k, k));
            }
            b = 1.0;
            for (i = 2; i <= matrix_height(D); i++) {
                b *= (matrix_get(D, i, i) - matrix_get(D, k, k));
            }
            cell = sign * sqrtl(a / (rho * b));
            matrix_set(v_hat, k, 1, cell);
        } else {
            long double num;
            long double denom;
            a = 1.0;
            for (i = 1; i <= k - 1; i++) {
                a *= (matrix_get(D, k, k) - matrix_get(lambda, i, 1));
            }
            b = 1.0;
            for (i = k; i <= n_deflated; i++) {
                b *= (matrix_get(lambda, i, 1) - matrix_get(D, k, k));
            }
            num = a * b;

            a = 1.0;
            for (i = 1; i <= k - 1; i++) {
                a *= (matrix_get(D, k, k) - matrix_get(D, i, i));
            }
            b = 1.0;
            for (i = k + 1; i <= n_deflated; i++) {
                b *= (matrix_get(D, i, i) - matrix_get(D, k, k));
            }

            denom = rho * a * b;
            
            cell = sign * sqrtl(num / denom);
            matrix_set(v_hat, k, 1, cell);
        }
    }

#ifdef DIVIDE_AND_CONQUER_DEBUG
    printf("v_hat:\n");
    matrix_print(v_hat);
#endif  /* DIVIDE_AND_CONQUER_DEBUG */
    
    jj = 1;
    for (i = 1; i <= n; i++) {
        double cell;
        unsigned kk;
        if (is_zero(matrix_get(v, i, 1))) {
            continue;
        }
        cell = matrix_get(lambda, jj, 1);
        matrix_set(eigenvalues, i, i, cell);
        tmp1 = matrix_diag(n_deflated, cell);
        tmp2 = matrix_minus(tmp1, D);
        matrix_free(tmp1);
        tmp1 = matrix_inverse(tmp2);
        matrix_free(tmp2);
        tmp2 = matrix_mul(tmp1, v_hat);
        matrix_free(tmp1);
        jj++;
        tmp1 = matrix_create(n, 1);
        matrix_set_def_val(tmp1, 0.0);
        kk = 1;
        for (j = 1; j <= n; j++) {
            if (is_zero(matrix_get(v, j, 1))) {
                matrix_set(tmp1, j, 1, 0.0);
            } else {
                cell = matrix_get(tmp2, kk, 1);
                matrix_set(tmp1, j, 1, cell);
                kk++;
            }
        }
        for (j = 1; j <= n; j++) {
            cell = matrix_get(tmp1, j, 1) / vector_norm(tmp1);
            matrix_set(eigenvectors, j, i, cell);
        }

        matrix_free(tmp1);
        matrix_free(tmp2);
    }

    matrix_free(v_hat);

    tmp1 = matrix_mul(G, eigenvectors);
    matrix_free(eigenvectors);
    eigenvectors = tmp1;

#ifdef DIVIDE_AND_CONQUER_DEBUG
    printf("eigenvectors:\n");
    matrix_print(eigenvectors);
    printf("eigenvalues:\n");
    matrix_print(eigenvalues);
#endif  /* DIVIDE_AND_CONQUER_DEBUG */    
    
    P1 = matrix_diag_permut(eigenvalues);
    tmp1 = matrix_mul(P1, eigenvalues);
    matrix_free(eigenvalues);
    tmp2 = matrix_transpose(P1);
    eigenvalues = matrix_mul(tmp1, tmp2);
    matrix_free(tmp1);
    matrix_free(tmp2);
    tmp1 = matrix_mul(eigenvectors, P1);
    matrix_free(eigenvectors);
    eigenvectors = tmp1;
    tmp1 = matrix_transpose(P);
    tmp2 = matrix_mul(tmp1, eigenvalues);
    matrix_free(eigenvalues);
    eigenvalues = matrix_mul(tmp2, P);
    matrix_free(tmp2);
    tmp2 = matrix_mul(tmp1, eigenvectors);
    matrix_free(eigenvectors);
    eigenvectors = matrix_mul(tmp2, P);
    matrix_free(tmp2);
    matrix_free(tmp1);
    (*L) = eigenvalues;


    (*Q) = matrix_from_two_blocks(Q1, Q2);
    tmp1 = matrix_mul((*Q), eigenvectors);
    matrix_free((*Q));
    (*Q) = tmp1;

#ifdef DIVIDE_AND_CONQUER_DEBUG
    printf("----divide-and-conquer-end----\n");
    printf("\n");
#endif  /* DIVIDE_AND_CONQUER_DEBUG */
    
    matrix_free(T1);
    matrix_free(Q1);
    matrix_free(L1);
    
    matrix_free(T2);
    matrix_free(Q2);
    matrix_free(L2);

    matrix_free(eigenvectors);
    matrix_free(v);
    matrix_free(v_prime);
    matrix_free(lambda);
    matrix_free(G);
    matrix_free(D);

    matrix_free(P);
    matrix_free(P1);
    
    return 0;

 err3:
    matrix_free(v);
 err2:
    matrix_free(D);
 err1:
    matrix_free(T2);
 err0:
    matrix_free(T1);
    return r;
}
