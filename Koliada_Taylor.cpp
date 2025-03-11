#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define PI M_PI

void sinx(float* terms, int n, float x)
{
    x = fmodf(x, 2.0f * (float)PI);
    terms[0] = x;
    for (int i = 1; i < n; i++) {
        terms[i] = (terms[i - 1] * -1 * x * x) / ((2 * i + 1) * (2 * i));
    }
}

void cosx(float* terms, int n, float x)
{
    x = fmodf(x, 2.0f * (float)PI);
    terms[0] = 1.0f;
    for (int i = 1; i < n; i++) {
        terms[i] = (terms[i - 1] * -1 * x * x) / ((2 * i - 1) * (2 * i));
    }
}

void ln(float* terms, int n, float x)
{
    terms[0] = x;
    for (int i = 1; i < n; i++) {
        terms[i] = (terms[i - 1] * -1 * x * i) / (i + 1);
    }
}

void expx(float* terms, int n, float x)
{
    terms[0] = 1.0f;
    for (int i = 1; i < n; i++) {
        terms[i] = (terms[i - 1] * x) / i;
    }
}

// Прямое суммирование
float direct_sum(float* terms, int n)
{
    float fx = 0.0f;
    for (int i = 0; i < n; i++) {
        fx += terms[i];
    }
    return fx;
}

// Обратное суммирование
float reverse_sum(float* terms, int n)
{
    float fx = 0.0f;
    for (int i = n - 1; i >= 0; i--) {
        fx += terms[i];
    }
    return fx;
}

// Попарное суммирование
float pair_sum(float* terms, int n)
{
    float fx = 0.0f;
    for (int i = 0; i < n / 2; i++) {
        fx += terms[2 * i];
        if (2 * i + 1 < n) {
            fx += terms[2 * i + 1];
        }
    }
    if (n % 2 != 0) {
        fx += terms[n - 1];
    }
    return fx;
}

void filling(float* terms, int n, int func, float x,
    float* direct_fx, float* reverse_fx,
    float* pair_fx, float* expected_res)
{
    if (func == 1) {
        sinx(terms, n, x);
        *expected_res = sin(x);
    }
    else if (func == 2) {
        cosx(terms, n, x);
        *expected_res = cos(x);
    }
    else if (func == 3) {
        ln(terms, n, x);
        *expected_res = log(1 + x);
    }
    else {
        expx(terms, n, x);
        *expected_res = exp(x);
    }

    *direct_fx = direct_sum(terms, n);
    *reverse_fx = reverse_sum(terms, n);
    *pair_fx = pair_sum(terms, n);
}

int main()
{
    int n, func;
    float x, direct_fx, reverse_fx, pair_fx, expected_res;

    printf("Select function:\n");
    printf("1. sinx\n");
    printf("2. cosx\n");
    printf("3. ln(1+x)\n");
    printf("4. e^x\n");
    scanf_s("%d", &func);

    printf("Enter number of elements\n");
    scanf_s("%d", &n);

    printf("Enter x for function\n");
    scanf_s("%f", &x);
    float* terms = (float*)malloc(sizeof(float) * n);
    filling(terms, n, func, x, &direct_fx, &reverse_fx, &pair_fx, &expected_res);

    printf("\nResults:\n");
    printf("Direct summation:  %f\n", direct_fx);
    printf("Reverse summation: %f\n", reverse_fx);
    printf("Pair summation:    %f\n", pair_fx);
    printf("Expected result:   %f\n", expected_res);
    printf("Max difference:    %f\n", fmaxf(fmaxf(fabs(expected_res - direct_fx), fabs(expected_res - reverse_fx)), fabs(expected_res - pair_fx)));

    free(terms);
    return 0;
}