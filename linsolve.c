#include "linsolve.h"
#include "error.h"
#include <stdlib.h>
#include <math.h>

inline static void swap(int* a, int* b)
{
	if (a == b)
	{
		return;
	}
	*a ^= *b;
	*b ^= *a;
	*a ^= *b;
}

int linsolve(int n, double* A, const double* b, int ldb, double* x)
{
	if (n <= 0)
	{
		return 0;
	}

	int* perm = (int*)malloc(sizeof(int) * n);
	CHK_MALLOC(perm);

	double* b_ = (double*)malloc(sizeof(double) * n);
	CHK_MALLOC(b_);

	for (int i = 0; i < n; i++)
	{
		perm[i] = i;
		b_[i] = b[ldb * i];
	}

#define A(i, j) A[perm[(i)] + (j) * n]
#define b(i) b_[perm[(i)]]
	for (int i = 0; i < n; i++)
	{
		int pivot = i;
		for (int k = i + 1; k < n; k++)
		{
			if (fabs(A(k, i)) > fabs(A(pivot, i)))
			{
				pivot = k;
			}
		}
		swap(perm + i, perm + pivot);

		double q = A(i, i);
		A(i, i) = 1;
		for (int j = i + 1; j < n; j++)
		{
			A(i, j) /= q;
		}
		b(i) /= q;

		for (int k = i + 1; k < n; k++)
		{
			double m = A(k, i);
			A(k, i) = 0.0;
			for (int j = i + 1; j < n; j++)
			{
				A(k, j) -= m * A(i, j);
			}
			b(k) -= m * b(i);
		}
	}

	for (int i = n - 2; i >= 0; i--)
	{
		for (int j = i + 1; j < n; j++)
		{
			b(i) -= A(i, j) * b(j);
		}
	}

	for (int i = 0; i < n; i++)
	{
		x[i] = b(i);
	}

	free(perm);
	return 0;
}