#include <stdio.h>
#include <math.h>

double psin(double x) {
    return x * (1 - 29593.0 / 207636.0 * x * x + 34911.0 / 7613320.0 * x*x*x*x - 479249.0 / 11511339840.0 * x*x*x*x*x*x ) / (1 + 1671.0 / 69212.0 * x * x + 97.0 / 351384.0 * x*x*x*x + 2623.0 / 1644477120.0 * x*x*x*x*x*x);
}

double tsin(double x){
    return x-x*x*x/6+x*x*x*x*x/120-x*x*x*x*x*x*x/5040;
}

void t3() {
    double pi = 3.1415926535897932;
    double x[5] = {pi / 2024, pi / 10, pi / 6, pi / 4, pi / 3};
    for(int i = 0; i < 5; i++) {
        if(fabs(x[i]) < pow(10, -8))
            printf("x = %f, psin1(x) = %.11E, tsin(x) = %.11E \n", x[i], x[i], tsin(x[i]));
        else if(fabs(x[i]) < (pi / 6))
            printf("x = %f, psin2(x) = %.11E & %.11E \n", x[i], psin(x[i]), tsin(x[i]));
        else
            printf("x = %f, psin3(x) = %.11E & %.11E \n", x[i], psin(x[i] / 3) * (3 - 4 * pow(sin(x[i] / 3), 2)), tsin(x[i]));
    }
}

void t2() {
    double a[7]={4042.045051380452, 0.000531415926535, -2759471.276702747, 0.0000557052996742895, 2755463.874010974, -34.64291531256604,-0.000031415926535};
    double sum[4];
    for (int i = 0; i < 7; ++i) {
        sum[0] += a[i], sum[1] += a[6 - i];
    }
    double b[7]={4042.045051380452, 0.000531415926535, -2759471.276702747, 0.0000557052996742895, 2755463.874010974, -34.64291531256604,-0.000031415926535};
    for (int i = 0; i < 7; ++i) {
        for(int j = i + 1; j < 7; ++j) {
            if (abs(b[i]) < abs(b[j])) {
                double temp = b[i];
                b[i] = b[j];
                b[j] = temp;
            }
        }
    }
    for(int i=0;i<7;i++) {
        printf("b[%d]=%.12E\n", i, b[i]);
    }
    for (int i = 0; i < 7; ++i) {
        sum[2] += b[i], sum[3] += b[6 - i];
    }
    for(int i=0;i<4;i++) {
        printf("sum[%d]=%.12E\n", i, sum[i]);
    }
}

void t1() {
    float x, a, b;
    for (int i = 0; i < 12; i++) {
        x = pow(4, -1 * i);
        a = (float)sqrt(x * x + (float)49.0) - (float) 7.0;
        b = x * x / ((float)sqrt(x * x + (float)49.0) + (float)7.0);
        printf("\\hline\n$4^{-%d}$ & %.12E & %.12E \\\\\n", i, a, b);
    }
}

int main() {
    t1();
    t2();
    t3();
}
