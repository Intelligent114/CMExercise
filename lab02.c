#include <stdio.h>
#include <math.h>

double f(double x) {
    return (3.0*x + 1.0) / (x*x - 2.0*x + 3.0);
}

double lagrange(double x, const double *X, const double *Y, int N) {
    double result = 0.0;
    for(int j = 0; j <= N; j++) {
        double L = 1.0;
        for(int m = 0; m <= N; m++) {
            if(m != j) {
                L *= (x - X[m]) / (X[j] - X[m]);
            }
        }
        result += Y[j] * L;
    }
    return result;
}

int main() {
    int Ns[4] = {4, 8, 16, 32};

    for(int k = 0; k < 4; k++) {
        int N = Ns[k];

        double xnode[200], ynode[200];

        for(int i = 0; i <= N; i++) {
            xnode[i] = -5.0 + 10.0 * i / N;
            ynode[i] = f(xnode[i]);
        }

        double maxErr = 0.0;
        for(int i = 0; i <= 500; i++) {
            double x = i / 50.0 - 5.0;
            double val = lagrange(x, xnode, ynode, N);
            double err = fabs(f(x) - val);
            if(err > maxErr) {
                maxErr = err;
            }
        }

        printf("\\hline\n$%d$ & $%.12e$", N, maxErr);

        for(int i = 0; i <= N; i++) {
            xnode[i] = -5.0 * cos((2.0 * i + 1.0) * M_PI / (2.0 * (N + 1)));
            ynode[i] = f(xnode[i]);
        }

        maxErr = 0.0;
        for(int i = 0; i <= 500; i++) {
            double x = i / 50.0 - 5.0;
            double val = lagrange(x, xnode, ynode, N);
            double err = fabs(f(x) - val);
            if(err > maxErr) {
                maxErr = err;
            }
        }
        printf(" & $%.12e$ \\\\\n", maxErr);
    }

    return 0;
}
