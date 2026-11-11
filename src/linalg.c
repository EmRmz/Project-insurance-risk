#include "linalg.h"

#include <math.h>

int cholesky(double *a, int k)
{
    for (int j = 0; j < k; j++) {
        double d = a[j * k + j];
        for (int p = 0; p < j; p++)
            d -= a[j * k + p] * a[j * k + p];
        if (d <= 0)
            return 0;
        a[j * k + j] = sqrt(d);
        for (int i = j + 1; i < k; i++) {
            double s = a[i * k + j];
            for (int p = 0; p < j; p++)
                s -= a[i * k + p] * a[j * k + p];
            a[i * k + j] = s / a[j * k + j];
        }
    }
    return 1;
}

void cholesky_solve(const double *l, int k, double *b)
{
    for (int i = 0; i < k; i++) {        /* L y = b */
        for (int p = 0; p < i; p++)
            b[i] -= l[i * k + p] * b[p];
        b[i] /= l[i * k + i];
    }
    for (int i = k - 1; i >= 0; i--) {   /* L' x = y */
        for (int p = i + 1; p < k; p++)
            b[i] -= l[p * k + i] * b[p];
        b[i] /= l[i * k + i];
    }
}

int gauss_solve(double *a, double *b, int k)
{
    for (int c = 0; c < k; c++) {
        int piv = c;
        for (int r = c + 1; r < k; r++)
            if (fabs(a[r * k + c]) > fabs(a[piv * k + c]))
                piv = r;
        if (fabs(a[piv * k + c]) < 1e-14)
            return 0;
        if (piv != c) {
            for (int j = 0; j < k; j++) {
                double t = a[c * k + j];
                a[c * k + j] = a[piv * k + j];
                a[piv * k + j] = t;
            }
            double t = b[c];
            b[c] = b[piv];
            b[piv] = t;
        }
        for (int r = c + 1; r < k; r++) {
            double f = a[r * k + c] / a[c * k + c];
            for (int j = c; j < k; j++)
                a[r * k + j] -= f * a[c * k + j];
            b[r] -= f * b[c];
        }
    }
    for (int r = k - 1; r >= 0; r--) {
        for (int j = r + 1; j < k; j++)
            b[r] -= a[r * k + j] * b[j];
        b[r] /= a[r * k + r];
    }
    return 1;
}
