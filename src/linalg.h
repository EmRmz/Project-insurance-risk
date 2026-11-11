/* Small dense linear algebra for k x k systems, row-major storage. */
#ifndef LINALG_H
#define LINALG_H

/* Cholesky factorization A = L L' in place (lower triangle). Returns 0 if A is not positive definite. */
int cholesky(double *a, int k);

/* Solve L L' x = b using the factor from cholesky(); b is overwritten with x. */
void cholesky_solve(const double *l, int k, double *b);

/* Solve A x = b by Gaussian elimination with partial pivoting; a and b are overwritten. */
int gauss_solve(double *a, double *b, int k);

#endif
