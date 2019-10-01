#include "deflation.h"
#include "matrix_algo.h"
#include "utils.h"
#include "macro.h"

#include <assert.h>
#include <stdlib.h>

/*
  If DEFLATION_DEBUG is defined,
  deflate will print
  some debug info to stdout.
 */
/* #define DEFLATION_DEBUG */

#ifdef DEFLATION_DEBUG
#include <stdio.h>
#endif  /* DEFLATION_DEBUG */

static matrix_type_t matrix_remove_ith(const matrix_type_t mat, unsigned k)
{
    unsigned i, j, l, m;
    unsigned n;

    assert(matrix_height(mat) == matrix_width(mat));

    n = matrix_height(mat);
    
    matrix_type_t res = matrix_create(n - 1, n - 1);
    if (res == NULL) {
        goto err0;
    }

    l = 1;

    for (i = 1; i <= n; i++) {
        if (i == k) {
            continue;
        }
        m = 1;
        for (j = 1; j <= n; j++) {
            long double cell;
            if (j == k) {
                continue;
            }

            cell = matrix_get(mat, i, j);
            matrix_set(res, l, m, cell);
            m++;
        }
        l++;
    }

    return res;

 err0:
    return NULL;
}

static void rotate(long double x1, long double x2, long double*c, long double*s)
{
    long double norm = sqrtl(x1 * x1 + x2 * x2);
    (*c) = x1 / norm;
    (*s) = -x2 / norm;
}

static void e(matrix_type_t eigenvectors, unsigned i)
{
    matrix_set(eigenvectors, i, i, 1.0);
}

int deflate(matrix_type_t*D, matrix_type_t*v, matrix_type_t*v_prime, matrix_type_t*eigenvalues, matrix_type_t*eigenvectors, unsigned*n_deflated, matrix_type_t*G, long double eps)
{
    int r;
    
    unsigned n;
    unsigned i, j;
    matrix_type_t tmp1, tmp2;
    unsigned*index_deflation;
    unsigned index_deflation_len;
    long double*v_prime_arr;
    unsigned v_prime_arr_len;
    unsigned ii, nn;

    assert(matrix_height((*D)) == matrix_width((*D)));
    assert((matrix_height((*v)) == matrix_height((*D))) && (matrix_width((*v)) == 1));

#ifdef DEFLATION_DEBUG
    printf("--------deflation-beg---------\n");
#endif  /* DEFLATION_DEBUG */

    n = matrix_height((*D));

    (*G) = matrix_diag(n, 1.0);
    if ((*G) == NULL) {
        r = -1;
        goto err0;
    }

    /*
      Compute plane rotation matrix G.
     */
    for (j = 1; j <= n - 1; j++) {
        if (are_equal(matrix_get((*D), j, j), matrix_get((*D), j + 1, j + 1))) {
            long double c, s;
            matrix_type_t G_tmp;
            rotate(matrix_get((*v), j, 1), matrix_get((*v), j + 1, 1), &c, &s);
            G_tmp = matrix_diag(n, 1.0);
            if (G_tmp == NULL) {
                r = -2;
                goto err1;
            }
            matrix_set(G_tmp, j, j, c);
            matrix_set(G_tmp, j, j + 1, s);
            matrix_set(G_tmp, j + 1, j, -s);
            matrix_set(G_tmp, j + 1, j + 1, c);
            tmp1 = matrix_mul(G_tmp, (*G));
            if (tmp1 == NULL) {
                matrix_free(G_tmp);
                r = -3;
                goto err1;
            }
            matrix_free(G_tmp);
            matrix_free((*G));
            (*G) = tmp1;
        }
    }

    /*
      Simplifying D and V.
     */
    (*n_deflated) = 0;
    index_deflation = (unsigned*) malloc(n * sizeof(unsigned));
    index_deflation_len = 0;

    tmp1 = matrix_transpose((*G));
    if (tmp1 == NULL) {
        r = -4;
        goto err1;
    }
    tmp2 = matrix_mul(tmp1, (*v));
    if (tmp2 == NULL) {
        matrix_free(tmp1);
        r = -5;
        goto err1;
    }
    matrix_free(tmp1);
    matrix_free((*v));
    (*v) = tmp2;

#ifdef DEFLATION_DEBUG
    printf("G:\n");
    matrix_print((*G));
    printf("v:\n");
    matrix_print((*v));
#endif  /* DEFLATION_DEBUG */

    v_prime_arr = (long double*) malloc(n * sizeof(long double));
    if (v_prime_arr == NULL) {
        r = -6;
        goto err2;
    }
    v_prime_arr_len = 0;

    for (i = 1; i <= n; i++) {
        if (is_zero(matrix_get((*v), i, 1))) {
            index_deflation[index_deflation_len++] = i;
        } else {
            (*n_deflated)++;
            v_prime_arr[v_prime_arr_len++] = matrix_get((*v), i, 1);
        }
    }

    (*v_prime) = matrix_create(1, v_prime_arr_len);
    if ((*v_prime) == NULL) {
        r = -7;
        goto err3;
    }
    for (i = 0; i < v_prime_arr_len; i++) {
        matrix_set((*v_prime), 1, i + 1, v_prime_arr[i]);
    }

    free(v_prime_arr);
    v_prime_arr = NULL;

#ifdef DEFLATION_DEBUG
    printf("n_deflated:\n");
    printf("%u\n", (*n_deflated));
    printf("v_prime:\n");
    matrix_print((*v_prime));
#endif  /* DEFLATION_DEBUG */

    ii = 1;
    nn = n;
    i = 1;

    matrix_type_t DD = matrix_copy((*D));
    if (DD == NULL) {
        r = -8;
        goto err4;
    }

#ifdef DEFLATION_DEBUG
    printf("D before deflation:\n");
    matrix_print((*D));
#endif  /* DEFLATION_DEBUG */

    /*
      Remove entries of D corresponding to deflated entries 
    */
    while (ii <= nn) {
        if (is_zero(matrix_get((*v), i, 1))) {
            tmp1 = matrix_remove_ith((*D), ii);
            if (tmp1 == NULL) {
                r = -9;
                goto err5;
            }
            matrix_free((*D));
            (*D) = tmp1;
            nn--;
        } else {
            ii++;
        }
        i++;
    }

#ifdef DEFLATION_DEBUG
    printf("D after deflation:\n");
    matrix_print((*D));
#endif  /* DEFLATION_DEBUG */    

    (*eigenvectors) = matrix_create(n, n);
    if ((*eigenvectors) == NULL) {
        r = -10;
        goto err5;
    }
    matrix_set_def_val((*eigenvectors), 0.0);
    (*eigenvalues) = matrix_create(n, n);
    if ((*eigenvalues) == NULL) {
        r = -11;
        goto err6;
    }
    matrix_set_def_val((*eigenvalues), 0.0);

    /*
      Compute eigevalues and eigenvectors for deflated cases.
    */
    ii = 1;
    for (i = 1; i <= n; i++) {
        if (is_zero(matrix_get((*v), i, 1))) {
            long double cell = matrix_get(DD, index_deflation[ii - 1], index_deflation[ii - 1]);
            matrix_set((*eigenvalues), i, i, cell);
            ii++;
            e((*eigenvectors), i);
        }
    }

    matrix_free(DD);
    free(index_deflation);

#ifdef DEFLATION_DEBUG
    printf("--------deflation-end---------\n");
#endif  /* DEFLATION_DEBUG */    

    return 0;

 err6:
    matrix_free((*eigenvectors));
 err5:
    matrix_free(DD);
 err4:
    matrix_free((*v_prime));
 err3:
    if (v_prime_arr != NULL) {
        free(v_prime_arr);
    }
 err2:
    matrix_free((*v));
 err1:
    matrix_free((*G));
 err0:
    return r;
}
