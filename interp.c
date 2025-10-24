#include "interp.h"
#include "linsolve.h"
#include "error.h"
#include <stdlib.h>
#include <math.h>

// Binary search of a cell in the time grid for the given time point 
inline static int findCell(int n, const double* ts, double t)
{
	if (n <= 2)
	{
		return 0;
	}

	int i = 0;
	if (t < ts[1])
	{
		i = 0;
	}
	else if (t >= ts[n - 2])
	{
		i = n - 2;
	}
	else
	{
		int a = 0;
		int b = n - 1;
		while (b - a > 1)
		{
			int c = (a + b) / 2;
			if (t < ts[c])
			{
				b = c;
			}
			else
			{
				a = c;
			}
		}
		i = a;
	}
	return i;
}

double lerp(int n, const double* x, const double* y, int ldy, double xq)
{
	if (n <= 0)
	{
		return 0.0;
	}
	if (n == 1)
	{
		return y[0];
	}

	int i = findCell(n, x, xq);

	double a = x[i];
	double b = x[i + 1];
	double A = y[ldy * i];
	double B = y[ldy * (i + 1)];
	double sigma = (xq - a) / (b - a);
	sigma = min(max(sigma, 0.0), 1.0);
	double yq = (1 - sigma) * A + sigma * B;
	return yq;
}

double lerpDer(int n, const double* x, const double* y, int ldy, double xq)
{
	if (n <= 1)
	{
		return 0.0;
	}

	int i = findCell(n, x, xq);

	double a = x[i];
	double b = x[i + 1];
	double A = y[ldy * i];
	double B = y[ldy * (i + 1)];
	double der = (B - A) / (b - a);
	return der;
}

double calcTrigPoly(int m, double omega, const double* a, double t)
{
	int k = 2 * m + 1;
	double p = a[0];
	double theta = omega * t;
	for (int j = 1; j <= m; j++)
	{
		double co = cos(j * theta);
		double si = sin(j * theta);
		p += a[2 * j - 1] * co + a[2 * j] * si;
	}
	return p;
}

double calcTrigPolyDer(int m, double omega, const double* a, double t)
{
	int k = 2 * m + 1;
	double p = 0.0;
	double theta = omega * t;
	for (int j = 1; j <= m; j++)
	{
		double co = cos(j * theta);
		double si = sin(j * theta);
		p += j * omega * (a[2 * j] * co - a[2 * j - 1] * si);
	}
	return p;
}

double calcSlidingTrigInterp(int n, int k, double omega, const double* a, const double* ts, double t)
{
	if (n <= 0)
	{
		return 0;
	}
	if (n == 1)
	{
		return a[0];
	}
	int i = findCell(n, ts, t);
	int m = (k - 1) / 2;
	i = min(max(i - m, 0), n - k);
	double p = calcTrigPoly(m, omega, a + k * i, t);
	return p;
}

double calcSlidingTrigInterpDer(int n, int k, double omega, const double* a, const double* ts, double t)
{
	if (n <= 0)
	{
		return 0;
	}
	if (n == 1)
	{
		return a[0];
	}
	int i = findCell(n, ts, t);
	int m = (k - 1) / 2;
	i = min(max(i - m, 0), n - k);
	double p = calcTrigPolyDer(m, omega, a + k * i, t);
	return p;
}

int findTrigPoly(int m, double omega, const double* t, const double* p, int ldp, double* a)
{
	int k = 2 * m + 1;
	double* A = malloc(sizeof(double) * k * k);
	CHK_MALLOC(A);

#define A(i, j) A[(i) + (j) * k]
	for (int i = 0; i < k; i++)
	{
		double theta = omega * t[i];
		A[i] = 1;
		for (int j = 1; j <= m; j++)
		{
			A(i, 2 * j - 1) = cos(j * theta);
			A(i, 2 * j) = sin(j * theta);
		}
	}

	linsolve(k, A, p, ldp, a);

	free(A);
	return 0;
}

int findSlidingTrigInterp(int n, int k, double omega, const double* t, const double* p, int ldp, double* a)
{
	if (n <= 0)
	{
		return 0;
	}
	k = min(k, n);
	int m = (k - 1) / 2;
	for (int i = 0; i <= n - k; i++)
	{
		if (findTrigPoly(m, omega, t + i, p + ldp * i, ldp, a + k * i))
		{
			return 1;
		}
	}
	return 0;
}