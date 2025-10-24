#pragma once

// Solves a system of linear equations using Gaussian elimination
// (ldb is the stride of b; ld = leading dimension, as in FORTRAN)
int linsolve(int n, double* A, const double* b, int ldb, double* x);