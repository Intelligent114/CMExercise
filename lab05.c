#include <stdio.h>
#include <math.h>

void gauss_eliminate(double mat[4][5]) {
    int n = 4;
    for (int i = 0; i < n; i++) {
        int max_row = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(mat[k][i]) > fabs(mat[max_row][i])) {
                max_row = k;
            }
        }
        for (int k = i; k <= n; k++) {
            double temp = mat[i][k];
            mat[i][k] = mat[max_row][k];
            mat[max_row][k] = temp;
        }
        double pivot = mat[i][i];
        if (pivot == 0) {
            printf("矩阵奇异，无法求解\n");
            return;
        }
        for (int j = i; j <= n; j++) {
            mat[i][j] /= pivot;
        }
        for (int k = i + 1; k < n; k++) {
            double factor = mat[k][i];
            for (int j = i; j <= n; j++) {
                mat[k][j] -= factor * mat[i][j];
            }
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        for (int k = i - 1; k >= 0; k--) {
            double factor = mat[k][i];
            mat[k][i] -= factor * mat[i][i];
            mat[k][n] -= factor * mat[i][n];
        }
    }
}

int main() {
    double t[] = {0, 1, 2, 3, 4, 5, 6, 7};
    double y_total[] = {14.0011, 14.0541, 14.1008, 14.1212, 14.1260, 14.1175, 14.1000, 14.0541};
    double y_birth[] = {1723, 1523, 1465, 1202, 1062, 956, 902, 954};
    int n = 8;

    double sum_t[7] = {0};
    for (int i = 0; i < n; i++) {
        double ti = t[i];
        double pow_ti = 1.0;
        for (int j = 0; j < 7; j++) {
            sum_t[j] += pow_ti;
            pow_ti *= ti;
        }
    }

    double sum_yt_total[4] = {0};
    for (int i = 0; i < n; i++) {
        double ti = t[i];
        double yi = y_total[i];
        sum_yt_total[0] += yi;
        sum_yt_total[1] += yi * ti;
        sum_yt_total[2] += yi * ti * ti;
        sum_yt_total[3] += yi * ti * ti * ti;
    }

    double aug_total[4][5] = {
            {sum_t[0], sum_t[1], sum_t[2], sum_t[3], sum_yt_total[0]},
            {sum_t[1], sum_t[2], sum_t[3], sum_t[4], sum_yt_total[1]},
            {sum_t[2], sum_t[3], sum_t[4], sum_t[5], sum_yt_total[2]},
            {sum_t[3], sum_t[4], sum_t[5], sum_t[6], sum_yt_total[3]}
    };
    gauss_eliminate(aug_total);
    double a0 = aug_total[0][4], a1 = aug_total[1][4], a2 = aug_total[2][4], a3 = aug_total[3][4];

    double sum_yt_birth[4] = {0};
    for (int i = 0; i < n; i++) {
        double ti = t[i];
        double yi = y_birth[i];
        sum_yt_birth[0] += yi;
        sum_yt_birth[1] += yi * ti;
        sum_yt_birth[2] += yi * ti * ti;
        sum_yt_birth[3] += yi * ti * ti * ti;
    }

    double aug_birth[4][5] = {
            {sum_t[0], sum_t[1], sum_t[2], sum_t[3], sum_yt_birth[0]},
            {sum_t[1], sum_t[2], sum_t[3], sum_t[4], sum_yt_birth[1]},
            {sum_t[2], sum_t[3], sum_t[4], sum_t[5], sum_yt_birth[2]},
            {sum_t[3], sum_t[4], sum_t[5], sum_t[6], sum_yt_birth[3]}
    };
    gauss_eliminate(aug_birth);
    double b0 = aug_birth[0][4], b1 = aug_birth[1][4], b2 = aug_birth[2][4], b3 = aug_birth[3][4];

    int years[] = {2025, 2030, 2035};
    printf("总人口预测（亿）:\n");
    for (int i = 0; i < 3; i++) {
        int t_pred = years[i] - 2017;
        double y = a0 + a1*t_pred + a2*t_pred*t_pred + a3*t_pred*t_pred*t_pred;
        printf("%d年: %.3f\n", years[i], y);
    }

    printf("\n出生人口预测（万）:\n");
    for (int i = 0; i < 3; i++) {
        int t_pred = years[i] - 2017;
        double y = b0 + b1*t_pred + b2*t_pred*t_pred + b3*t_pred*t_pred*t_pred;
        printf("%d年: %.1f\n", years[i], y);
    }

    return 0;
}