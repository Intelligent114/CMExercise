#include <stdio.h>
#include <math.h>

double f(double x) {
    return atan(x) + 0.2 * x * sin(x / 2) + 0.601958;
}

double df(double x) {
    return 1.0 / (1 + x*x) + 0.2 * sin(x/2) + 0.1 * x * cos(x/2);
}

double d2f(double x) {
    return -2*x / pow(1 + x*x, 2) + 0.2 * cos(x/2) - 0.05 * x * sin(x/2);
}

int newton(double x0, double *result, int *steps) {
    double x = x0;
    double fx = f(x);
    *steps = 0;
    while (fabs(fx) >= 1e-8 && (*steps) <= 1e4) {
        double dfx = df(x);
        if (dfx == 0) return -1;
        x -= fx / dfx;
        fx = f(x);
        (*steps)++;
    }
    if (fabs(fx) < 1e-8) {
        *result = x;
        return 0;
    } else {
        return -1;
    }
}


int householder(double x0, double *result, int *steps) {
    double x = x0;
    double fx = f(x);
    *steps = 0;
    while (fabs(fx) >= 1e-8 && (*steps) <= 1e4) {
        double dfx = df(x);
        double d2fx = d2f(x);
        double denom = 2 * dfx * dfx - fx * d2fx;
        if (denom == 0) return -1;
        x -= (2 * fx * dfx) / denom;
        fx = f(x);
        (*steps)++;
    }
    if (fabs(fx) < 1e-8) {
        *result = x;
        return 0;
    } else {
        return -1;
    }
}

int main() {
    double x0_list[] = {-75, -60, -50, -40, -30, -20, -10, -5, 0, 6, 15, 25, 35, 45, 50, 60, 75};
    int n = sizeof(x0_list) / sizeof(double);

    printf("Newton Method:\n");
    printf("x0\tSteps\tSolution\n");
    for (int i = 0; i < n; i++) {
        double x0 = x0_list[i], result;
        int steps;
        if (newton(x0, &result, &steps) == 0) {
            printf("%.0f\t%d\t%.6f\n", x0, steps, result);
        } else {
            printf("%.0f\tFailed\n", x0);
        }
    }

    printf("\nHouseholder Method:\n");
    printf("x0\tSteps\tSolution\n");
    for (int i = 0; i < n; i++) {
        double x0 = x0_list[i], result;
        int steps;
        if (householder(x0, &result, &steps) == 0) {
            printf("%.0f\t%d\t%.6f\n", x0, steps, result);
        } else {
            printf("%.0f\tFailed\n", x0);
        }
    }

    return 0;
}
