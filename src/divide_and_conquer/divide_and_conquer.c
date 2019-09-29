#include "divide_and_conquer.h"

#include "macro.h"
#include "matrix_algo.h"
#include "deflation.h"
#include "secular_equation.h"

#include <stddef.h>
#include <assert.h>
#include <math.h>

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
            matrix_set(res, i - b, j - b, cell);
        }
    }

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
    matrix_type_t ind;
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
        cell = matrix_get(Q2, 1, i);
        matrix_set(v, i, 1, cell);
    }

    ind = matrix_diag_permut(D);
    tmp1 = matrix_mul(ind, D);
    tmp2 = matrix_transpose(ind);
    matrix_free(D);
    D = matrix_mul(tmp1, tmp2);
    matrix_free(tmp1);
    matrix_free(tmp2);

    tmp1 = matrix_mul(ind, v);
    matrix_free(v);
    v = tmp1;

    r = deflate(D, v, &v_prime, &eigenvalues, &eigenvectors, &n_deflated, &G, eps);
    if (r) {
        r = -9;
        goto err2;
    }

    lambda = matrix_create(n_deflated, 1);
    if (lambda == NULL) {
        r = -10;
        goto err3;
    }

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
        cell = solve_secular_equation(rho, D, matrix_transpose(v_prime), i, lambda_init, n_deflated, eps);
        matrix_set(lambda, i, 1, cell);
    }

    v_hat = matrix_create(n_deflated, 1);
    for (k = 1; k <= n_deflated; k++) {
        double cell;
        if (k == 1) {
            cell = 0.0; /* TODO */
            matrix_set(v_hat, k, 1, cell);
        } else {
            cell = 0.0; /* TODO */
            matrix_set(v_hat, k, 1, cell);
        }
    }
    jj = 1;
    for (i = 1; i <= n; i++) {
        double cell;
        if (is_zero(matrix_get(v, i, 1))) {
            continue;
        }
        cell = matrix_get(lambda, jj, 1);
        matrix_set(eigenvalues, i, i, cell);
        /* TODO. */
        jj++;
    }

    tmp1 = matrix_mul(G, eigenvectors);
    matrix_free(eigenvectors);
    eigenvectors = tmp1;
    

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
