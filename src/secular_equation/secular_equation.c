#include "secular_equation.h"
#include "macro.h"

#include <math.h>

static long double calculate_psi1(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i)
{
    unsigned k;
    
    long double res = 0.0;

    for (k = 0; k <= i; k++) {
        res += ((matrix_get(u, k, 0) * matrix_get(u, k, 0)) / (matrix_get(D, k, k) - lambda));
    }

    res *= ro;

    return res;
}

static long double calculate_psi1_der(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i)
{
    unsigned k;

    long double res = 0.0;

    for (k = 0; k <= i; k++) {
        res += ((matrix_get(u, k, 0) * matrix_get(u, k, 0)) / ((matrix_get(D, k, k) - lambda) * (matrix_get(D, k, k) - lambda)));
    }

    res *= ro;

    return res;
}

static long double calculate_psi2(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i)
{
    unsigned k;
    
    long double res = 0.0;

    for (k = i + 1; k < matrix_height(u); k++) {
        res += ((matrix_get(u, k, 0) * matrix_get(u, k, 0)) / (matrix_get(D, k, k) - lambda));
    }

    res *= ro;

    return res;
}

static long double calculate_psi2_der(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i)
{
    unsigned k;

    long double res = 0.0;

    for (k = i + 1; k < matrix_height(u); k++) {
        res += ((matrix_get(u, k, 0) * matrix_get(u, k, 0)) / ((matrix_get(D, k, k) - lambda) * (matrix_get(D, k, k) - lambda)));
    }

    res *= ro;

    return res;
}

static long double calculate_c1(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i)
{
    long double psi1_der = calculate_psi1_der(D, u, lambda, ro, i);
    long double res = psi1_der * ((matrix_get(D, i, i) - lambda) * (matrix_get(D, i, i) - lambda));
    return res;
}

static long double calculate_c1_circum(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i)
{
    long double psi1_der = calculate_psi1_der(D, u, lambda, ro, i);
    long double psi1 = calculate_psi1(D, u, lambda, ro, i);
    long double res = psi1 - psi1_der * (matrix_get(D, i, i) - lambda);
    return res;
}

static long double calculate_c2(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i)
{
    long double psi2_der = calculate_psi2_der(D, u, lambda, ro, i);
    long double res = psi2_der * ((matrix_get(D, i + 1, i + 1) - lambda) * (matrix_get(D, i + 1, i + 1) - lambda));
    return res;
}

static long double calculate_c2_circum(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i)
{
    long double psi2_der = calculate_psi2_der(D, u, lambda, ro, i);
    long double psi2 = calculate_psi2(D, u, lambda, ro, i);
    long double res = psi2 - psi2_der * (matrix_get(D, i + 1, i + 1) - lambda);
    return res;
}

static long double calculate_c3(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i) {
    long double c1_circum = calculate_c1_circum(D, u, lambda, ro, i);
    long double c2_circum = calculate_c2_circum(D, u, lambda, ro, i);    
    return (1 + c1_circum + c2_circum);
}

#include <stdio.h>

static int should_stop_iterations(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, long double eps)
{
    unsigned k;
    
    long double res = 0.0;

    for (k = 0; k < matrix_height(u); k++) {
        res += ((matrix_get(u, k, 0) * matrix_get(u, k, 0)) / (matrix_get(D, k, k) - lambda));
    }

    res *= ro;

    res += 1.0;

    return is_zero(res);
}

#include <stdio.h>
#include <stdlib.h>

static long double solve_secular_equation_d0_plus_inf_iter(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, long double eps)
{
    unsigned i = 0;
    
    long double c1 = calculate_c1(D, u, lambda, ro, i);
    long double c2 = calculate_c2(D, u, lambda, ro, i);

    long double lambda_1 = matrix_get(D, i, i) + (c2 / c1);

    if (lambda_1 > matrix_get(D, i, i)) {
        if (should_stop_iterations(D, u, lambda_1, ro, eps)) {
            return lambda_1;
        }
        return solve_secular_equation_d0_plus_inf_iter(D, u, lambda_1, ro, eps);
    }

    printf("lel\n");
    exit(1);
}

long double solve_secular_equation_d0_plus_inf(const matrix_type_t D, const matrix_type_t u, long double ro, long double eps)
{
    long double lambda = matrix_get(D, 0, 0);

    while (1) {
        if (should_stop_iterations(D, u, lambda, ro, eps)) {
            return lambda;
        }

        lambda += eps;
    }
    
    /* return solve_secular_equation_d0_plus_inf_iter(D, u, lambda_start, ro, eps); */
}

static long double solve_secular_equation_common_iter(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, unsigned i, long double eps)
{
    long double c1 = calculate_c1(D, u, lambda, ro, i);
    long double c2 = calculate_c2(D, u, lambda, ro, i);
    long double c3 = calculate_c3(D, u, lambda, ro, i);
    long double d_i = matrix_get(D, i, i);
    long double d_i_1 = matrix_get(D, i + 1, i + 1);

    long double a = c3;
    long double b = -(c1 + c2 + c3 * (d_i + d_i_1));
    long double c = c1 * d_i_1 + c2 * d_i + c3 * d_i * d_i_1;

    long double disc = b * b - 4 * a * c;

    long double lambda_1 = (-b + sqrt(disc)) / (2 * a);
    long double lambda_2 = (-b - sqrt(disc)) / (2 * a);

    if ((lambda_1 > d_i_1) && (lambda_1 < d_i)) {
        if (should_stop_iterations(D, u, lambda_1, ro, eps)) {
            return lambda_1;
        }
        return solve_secular_equation_common_iter(D, u, lambda_1, ro, i, eps);
    } else if ((lambda_2 > d_i_1) && (lambda_2 < d_i)) {
        if (should_stop_iterations(D, u, lambda_2, ro, eps)) {
            return lambda_2;
        }
        return solve_secular_equation_common_iter(D, u, lambda_2, ro, i, eps);
    }

    printf("kek\n");
    exit(1);
}

long double solve_secular_equation_common(const matrix_type_t D, const matrix_type_t u, long double ro, unsigned i, long double eps)
{
    long double lambda_start = lambda_start = (matrix_get(D, i, i) + matrix_get(D, i + 1, i + 1)) / 2.0;
    return solve_secular_equation_common_iter(D, u, lambda_start, ro, i, eps);
}

static long double solve_secular_equation_minus_inf_dn_iter(const matrix_type_t D, const matrix_type_t u, long double lambda, long double ro, long double eps)
{
    unsigned i = matrix_height(D) - 1;
    
    long double c1 = calculate_c1(D, u, lambda, ro, i);
    long double c2 = calculate_c2(D, u, lambda, ro, i);

    long double lambda_1 = matrix_get(D, i, i) + (c2 / c1);

    if (lambda_1 < matrix_get(D, i, i)) {
        if (should_stop_iterations(D, u, lambda, ro, eps)) {
            return lambda_1;
        }
        return solve_secular_equation_minus_inf_dn_iter(D, u, lambda, ro, eps);
    }

    printf("plel\n");
    exit(1);
}

long double solve_secular_equation_minus_inf_dn(const matrix_type_t D, const matrix_type_t u, long double ro, long double eps)
{
    long double lambda = matrix_get(D, matrix_height(D) - 1, matrix_height(D) - 1);

    while (1) {
        if (should_stop_iterations(D, u, lambda, ro, eps)) {
            return lambda;
        }

        lambda -= eps;
    }
}
