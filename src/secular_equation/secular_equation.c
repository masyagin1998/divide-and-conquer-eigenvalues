#include "secular_equation.h"

#include <math.h>

static long double calculate_c1(matrix_type_t D, matrix_type_t v, unsigned i, long double lambda_init)
{
    unsigned k;

    long double res = 0.0;

    for (k = 1; k <= i; k++) {
        res += ((matrix_get(v, 1, k) * matrix_get(v, 1, k)) / ((matrix_get(D, k, k) - lambda_init) * (matrix_get(D, k, k) - lambda_init)));
    }

    res *= ((matrix_get(D, i, i) - lambda_init) * (matrix_get(D, i, i) - lambda_init));

    return res;
}

static long double calculate_c1_hat(matrix_type_t D, matrix_type_t v, unsigned i, long double lambda_init)
{
    unsigned k;

    long double res = 0.0;

    return res;
}

static long double calculate_c2(matrix_type_t D, matrix_type_t v, unsigned i, long double lambda_init)
{
    unsigned k;
    unsigned n = matrix_height(D);

    long double res = 0.0;

    for (k = i + 1; k <= n; k++) {
        res += ((matrix_get(v, 1, k) * matrix_get(v, 1, k)) / ((matrix_get(D, k, k) - lambda_init) * (matrix_get(D, k, k) - lambda_init)));
    }

    res *= ((matrix_get(D, i + 1, i + 1) - lambda_init) * (matrix_get(D, i + 1, i + 1) - lambda_init));

    return res;
}

static long double calculate_c2_hat(matrix_type_t D, matrix_type_t v, unsigned i, long double lambda_init)
{
    unsigned k;

    long double res = 0.0;

    return res;
}

static long double calculate_c3(long double c1_hat, long double c2_hat, long double rho)
{
    return 1 / rho + c1_hat + c2_hat;
}

static long double calculate_func(matrix_type_t D, matrix_type_t v, unsigned i, long double lambda_init)
{
    unsigned k;

    long double res = 0.0;

    return res;
}

long double solve_secular_equation(long double rho, matrix_type_t D, matrix_type_t v, unsigned i, long double lambda_init, unsigned n, long double eps)
{
    long double error = 1.0;
    unsigned iter = 1;

    while (error > eps) {
        long double c1 = calculate_c1(D, v, i, lambda_init);
        long double c1_hat = calculate_c1_hat(D, v, i, lambda_init);
        if (i == n) {
            lambda_init = (c1 * rho + c1_hat * rho * matrix_get(D, n, n) + matrix_get(D, n, n)) / (1 + rho * c1_hat);
            error = fabsl(1 / rho + calculate_func(D, v, i, lambda_init));
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

            error = fabsl(1 / rho /* + todo. */ );
        }
        iter++;
    }

    return lambda_init;
}
