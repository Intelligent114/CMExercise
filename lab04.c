
//代码运行时间较长，请耐心等待

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void compute_a_matrix(double *a, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i * n + j] = 1.0 / (i + j + 1);
        }
    }
}

void compute_diag(double *diag, int n) {
    for (int i = 0; i < n; i++) {
        diag[i] = 1.0 / (2 * i + 1) + 2.25;
    }
}

void compute_b_vector(double *b, double *a, int n) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += a[i * n + j];
        }
        b[i] = sum + 2.25;
    }
}

double jacobi(double *a, double *diag, double *b, int n, int max_iter, double epsilon, int *iterations) {
    double *x_old = (double *)calloc(n, sizeof(double));
    double *x_new = (double *)calloc(n, sizeof(double));
    double error = 0.0;
    *iterations = 0;

    if (!x_old || !x_new) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    while (*iterations < max_iter) {
        (*iterations)++;
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++) {
                if (j != i) sum += a[i * n + j] * x_old[j];
            }
            x_new[i] = (b[i] - sum) / diag[i];
        }

        double diff = 0.0;
        for (int i = 0; i < n; i++) {
            diff += fabs(x_new[i] - x_old[i]);
        }
        if (diff < epsilon) break;

        memcpy(x_old, x_new, n * sizeof(double));
    }

    for (int i = 0; i < n; i++) {
        error += fabs(x_new[i] - 1.0);
    }

    free(x_old);
    free(x_new);
    return error;
}

double gauss_seidel(double *a, double *diag, double *b, int n, int max_iter, double epsilon, int *iterations) {
    double *x = (double *)calloc(n, sizeof(double));
    double *x_old = (double *)malloc(n * sizeof(double));
    double error = 0.0;
    *iterations = 0;

    if (!x || !x_old) {
        fprintf(stderr, "Memory allocation failed.\n");
        exit(1);
    }

    while (*iterations < max_iter) {
        (*iterations)++;
        memcpy(x_old, x, n * sizeof(double));

        for (int i = 0; i < n; i++) {
            double sum_lower = 0.0;
            for (int j = 0; j < i; j++) {
                sum_lower += a[i * n + j] * x[j];
            }
            double sum_upper = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum_upper += a[i * n + j] * x_old[j];
            }
            x[i] = (b[i] - sum_lower - sum_upper) / diag[i];
        }

        double diff = 0.0;
        for (int i = 0; i < n; i++) {
            diff += fabs(x[i] - x_old[i]);
        }
        if (diff < epsilon) break;
    }

    for (int i = 0; i < n; i++) {
        error += fabs(x[i] - 1.0);
    }

    free(x);
    free(x_old);
    return error;
}

int main() {
    int n_values[] = {10, 30, 100, 500, 1500, 5000};
    int num_n = sizeof(n_values) / sizeof(n_values[0]);
    double epsilon = 1e-5;
    int max_iter = 500000;

    for (int i = 0; i < num_n; i++) {
        int n = n_values[i];
        printf("n = %d\n", n);

        double *a = (double *)malloc(n * n * sizeof(double));
        double *diag = (double *)malloc(n * sizeof(double));
        double *b = (double *)malloc(n * sizeof(double));

        if (!a || !diag || !b) {
            fprintf(stderr, "Memory allocation failed for n = %d.\n", n);
            continue;
        }

        compute_a_matrix(a, n);
        compute_diag(diag, n);
        compute_b_vector(b, a, n);

        int iter_jacobi = 0, iter_gs = 0;
        double error_jacobi = jacobi(a, diag, b, n, max_iter, epsilon, &iter_jacobi);
        double error_gs = gauss_seidel(a, diag, b, n, max_iter, epsilon, &iter_gs);

        printf("Jacobi Error: %.2e, Iterations: %d\n", error_jacobi, iter_jacobi);
        printf("Gauss-Seidel Error: %.2e, Iterations: %d\n\n", error_gs, iter_gs);

        free(a);
        free(diag);
        free(b);
    }

    return 0;
}