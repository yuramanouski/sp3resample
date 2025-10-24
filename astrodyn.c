#include "astrodyn.h"
#include "interp.h"
#include "error.h"
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

static const double omegaEarth = M_PI / (12 * 3600);
static const double muEarth = 3.986004418e14;

inline static rotz(double theta, const double* r, double* rOut)
{
	double co = cos(theta);
	double si = sin(theta);
	rOut[0] = co * r[0] - si * r[1];
	rOut[1] = si * r[0] + co * r[1];
	rOut[2] = r[2];
}

static inline double dot(const double* u, const double* v)
{
	double d = u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
	return d;
}

static inline double len(const double* u)
{
	double r = sqrt(dot(u, u));
	return r;
}

ecef2eci(double t, const double* rEcef, double* rEci, const double* vEcef, double* vEci)
{
	double theta = omegaEarth * t;
	rotz(theta, rEcef, rEci);
	if (vEci != NULL)
	{
		rotz(theta, vEcef, vEci);
		vEci[0] -= omegaEarth * rEci[1];
		vEci[1] += omegaEarth * rEci[0];
	}
}

eci2ecef(double t, const double* rEci, double* rEcef, const double* vEci, double* vEcef)
{
	double theta = -omegaEarth * t;
	rotz(theta, rEci, rEcef);
	if (vEci != NULL)
	{
		rotz(theta, vEci, vEcef);
		vEcef[0] += omegaEarth * rEcef[1];
		vEcef[1] -= omegaEarth * rEcef[0];
	}
}

ecef2eciN(int n, const double* ts, const double* rEcef, int ldrEcef, double* rEci, int ldrEci,
	const double* vEcef, int ldvEcef, double* vEci, int ldvEci)
{
	for (int i = 0; i < n; i++)
	{
		double t = ts[i];
		ecef2eci(t, rEcef + ldrEcef * i, rEci + ldrEci * i,
			vEcef == NULL ? NULL : vEcef + ldvEcef * i,
			vEci == NULL ? NULL : vEci + ldvEci * i);
	}
}

eci2ecefN(int n, const double* ts, const double* rEci, int ldrEci, double* rEcef, int ldrEcef,
	const double* vEci, int ldvEci, double* vEcef, int ldvEcef)
{
	for (int i = 0; i < n; i++)
	{
		double t = ts[i];
		eci2ecef(t, rEci + ldrEci * i, rEcef + ldrEcef * i,
			vEci == NULL ? NULL : vEci + ldvEci * i,
			vEcef == NULL ? NULL : vEcef + ldvEcef * i);
	}
}

double calcAngFreq(double a)
{
	double omega = sqrt(muEarth / (a * a * a));
	return omega;
}

double calcMajSemiaxisFromTraj(int n, const double* r3dKm)
{
	double a = 0.0;
	for (int i = 0; i < n; i++)
	{
		double r = 1e3 * len(r3dKm + 3 * i);
		a = max(a, r);
	}
	return a;
}

double calcAngFreqFromTraj(int n, const double* r3dKm)
{
	double a = calcMajSemiaxisFromTraj(n, r3dKm);
	double omega = calcAngFreq(a);
	return omega;
}

int resampleSatTrajEci(int nEpochs, const double* ts, double* rEci, int ldrEci,
	int nEpochsOut, const double* tsOut, double* rEciOut, int ldrEciOut, double* vEciOut, int ldvEciOut)
{
	double omega = calcAngFreqFromTraj(nEpochs, rEci);

	const int k = 9;
	int nInterpCoefs = k * (nEpochs - k + 1);
	double* a = (double*)malloc(sizeof(double) * nInterpCoefs);
	CHK_MALLOC(a);

	for (int j = 0; j < 3; j++)
	{
		findSlidingTrigInterp(nEpochs, k, omega, ts, rEci + j, ldrEci, a);

		for (int iOut = 0; iOut < nEpochsOut; iOut++)
		{
			double tOut = tsOut[iOut];
			if (rEciOut != NULL)
			{
				double coordinate = calcSlidingTrigInterp(nEpochs, k, omega, a, ts, tOut);
				rEciOut[ldrEciOut * iOut + j] = coordinate;
			}
			if (vEciOut != NULL)
			{
				double velocityComponent = calcSlidingTrigInterpDer(nEpochs, k, omega, a, ts, tOut);
				velocityComponent *= 1e4;
				vEciOut[ldvEciOut * iOut + j] = velocityComponent;
			}
		}
	}

	return 0;
}

int resampleSatTrajEcef(int nEpochs, const double* ts, double* rEcef, int ldrEcef,
	int nEpochsOut, const double* tsOut, double* rEcefOut, int ldrEcefOut, double* vEcefOut, int ldvEcefOut)
{
	double* rEci = (double*)malloc(sizeof(double) * 3 * nEpochs);
	CHK_MALLOC(rEci);

	double* rEciOut = (double*)malloc(sizeof(double) * 3 * nEpochsOut);
	CHK_MALLOC(rEciOut);

	double* vEciOut = NULL;
	
	if (vEcefOut != NULL)
	{
		vEciOut = (double*)malloc(sizeof(double) * 3 * nEpochsOut);
		CHK_MALLOC(vEciOut);
	}

	ecef2eciN(nEpochs, ts, rEcef, ldrEcef, rEci, 3, NULL, 0, NULL, 0);
	resampleSatTrajEci(nEpochs, ts, rEci, 3, nEpochsOut, tsOut, rEciOut, 3, vEciOut, 3);
	eci2ecefN(nEpochsOut, tsOut, rEciOut, 3, rEcefOut, ldrEcefOut, vEciOut, 3, vEcefOut, ldvEcefOut);

	free(rEci);
	free(rEciOut);
	free(vEciOut);
	return 0;
}

int resampleSatClock(int nEpochs, const double* ts, double* clock, int ldclock,
	int nEpochsOut, const double* tsOut, double* clockOut, int ldclockOut, double* clockRateOut, int ldclockRateOut)
{
	for (int iOut = 0; iOut < nEpochsOut; iOut++)
	{
		double tOut = tsOut[iOut];
		if (clockOut != NULL)
		{
			double c = lerp(nEpochs, ts, clock, ldclock, tOut);
			clockOut[ldclockOut * iOut] = c;
		}
		if (clockRateOut != NULL)
		{
			double cr = lerpDer(nEpochs, ts, clock, ldclock, tOut);
			cr *= 1e4;
			clockRateOut[ldclockRateOut * iOut] = cr;
		}
	}
	return 0;
}