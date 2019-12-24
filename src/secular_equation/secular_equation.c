#include "secular_equation.h"

#include <math.h>

/*
  If SECULAR_EQUATION_DEBUG is defined,
  secular_equation will print
  some debug info to stdout.
 */
/* #define SECULAR_EQUATION_DEBUG */

#ifdef SECULAR_EQUATION_DEBUG
#include <stdio.h>
#include "utils.h"
#endif  /* SECULAR_EQUATION_DEBUG */

static long double calculate_c1(const matrix_type_t D, const matrix_type_t v, unsigned i, long double lambda_init)
{
    unsigned k;

    long double res = 0.0;

    for (k = 1; k <= i; k++) {
        res += ((matrix_get(v, k, 1) * matrix_get(v, k, 1)) / ((matrix_get(D, k, k) - lambda_init) * (matrix_get(D, k, k) - lambda_init)));
    }

    res *= ((matrix_get(D, i, i) - lambda_init) * (matrix_get(D, i, i) - lambda_init));

    return res;
}

static long double calculate_c1_hat(const matrix_type_t D, const matrix_type_t v, unsigned i, long double lambda_init)
{
    unsigned k;

    long double res = 0.0;

    for (k = 1; k <= i; k++) {
        res += ((matrix_get(v, k, 1) * matrix_get(v, k, 1)) * (matrix_get(D, k, k) - matrix_get(D, i, i)) / ((matrix_get(D, k, k) - lambda_init) * (matrix_get(D, k, k) - lambda_init)));
    }

    return res;
}

static long double calculate_c2(const matrix_type_t D, const matrix_type_t v, unsigned i, long double lambda_init)
{
    unsigned k;
    unsigned n = matrix_height(D);

    long double res = 0.0;

    for (k = i + 1; k <= n; k++) {
        res += ((matrix_get(v, k, 1) * matrix_get(v, k, 1)) / ((matrix_get(D, k, k) - lambda_init) * (matrix_get(D, k, k) - lambda_init)));
    }

    res *= ((matrix_get(D, i + 1, i + 1) - lambda_init) * (matrix_get(D, i + 1, i + 1) - lambda_init));

    return res;
}

static long double calculate_c2_hat(const matrix_type_t D, const matrix_type_t v, unsigned i, long double lambda_init)
{
    unsigned k;
    unsigned n = matrix_height(D);

    long double res = 0.0;

    for (k = i + 1; k <= n; k++) {
        res += ((matrix_get(v, k, 1) * matrix_get(v, k, 1)) * (matrix_get(D, k, k) - matrix_get(D, i + 1, i + 1)) / ((matrix_get(D, k, k) - lambda_init) * (matrix_get(D, k, k) - lambda_init)));
    }

    return res;
}

static long double calculate_c3(long double c1_hat, long double c2_hat, long double rho)
{
    return 1 / rho + c1_hat + c2_hat;
}

static long double calculate_func(const matrix_type_t D, const matrix_type_t v, long double lambda_init)
{
    unsigned k;
    unsigned n = matrix_height(D);

    long double res = 0.0;

    for (k = 1; k <= n; k++) {
        res += ((matrix_get(v, k, 1) * matrix_get(v, k, 1)) / (matrix_get(D, k, k) - lambda_init));
    }

    return res;
}

long double solve_secular_equation(long double rho, const matrix_type_t D, const matrix_type_t v, unsigned i, long double lambda_init, unsigned n, long double eps)
{
    long double error = 1.0;
    unsigned iter = 1;

#ifdef SECULAR_EQUATION_DEBUG
    printf("-----secular-equation-beg-----\n");
#endif  /* SECULAR_EQUATION_DEBUG */

#ifdef SECULAR_EQUATION_DEBUG
    printf("rho:\n");
    printf("%Lf\n", rho);
    printf("D:\n");
    matrix_print(D, eps);
    printf("v:\n");
    matrix_print(v, eps);
    printf("i:\n");
    printf("%u\n", i);
    printf("lambda init:\n");
    printf("%Lf\n", lambda_init);
#endif  /* SECULAR_EQUATION_DEBUG */

    while (error > eps) {
        long double c1 = calculate_c1(D, v, i, lambda_init);
        long double c1_hat = calculate_c1_hat(D, v, i, lambda_init);
        if (i == n) {
            if ((1 + rho * c1_hat) != 0.0) {
                lambda_init = (c1 * rho + c1_hat * rho * matrix_get(D, n, n) + matrix_get(D, n, n)) / (1 + rho * c1_hat);
            } else {
                lambda_init = (c1 * rho + c1_hat * rho * matrix_get(D, n, n) + matrix_get(D, n, n)) / (eps / 1000.0);
            }
        } else {
            long double c2 = calculate_c2(D, v, i, lambda_init);
            long double c2_hat = calculate_c2_hat(D, v, i, lambda_init);
            long double c3 = calculate_c3(c1_hat, c2_hat, rho);
            long double p = -(c1 + c2 + c3 * matrix_get(D, i + 1, i + 1) + c3 * matrix_get(D, i, i)) / c3;
            long double q = (c1 * matrix_get(D, i + 1, i + 1) + c2 * matrix_get(D, i, i) + c3 * matrix_get(D, i, i) * matrix_get(D, i + 1, i + 1)) / c3;
            long double root1 = -p / 2.0 + sqrtl(p * p / 4 - q);
            long double root2 = -p / 2.0 - sqrtl(p * p / 4 - q);
            if ((matrix_get(D, i, i) <= root1) && (root1 <= matrix_get(D, i + 1, i + 1))) {
                lambda_init = root1;
            } else {
                lambda_init = root2;
            }
        }
        error = fabsl(1 / rho + calculate_func(D, v, lambda_init));        
        iter++;
        if (iter > 10) {
            break;
        }        
    }

#ifdef SECULAR_EQUATION_DEBUG
    printf("lambda:\n");
    printf("%Lf\n", lambda_init);
#endif  /* SECULAR_EQUATION_DEBUG */

#ifdef SECULAR_EQUATION_DEBUG
    printf("-----secular-equation-end-----\n");
#endif  /* SECULAR_EQUATION_DEBUG */

    return lambda_init;
}
