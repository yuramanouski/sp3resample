#pragma once
#include "interp.h"

// Converts a position and velocity from ECEF to ECI reference frame for a given time
ecef2eci(double t, const double* rEcef, double* rEci, const double* vEcef, double* vEci);

// Converts a position and velocity from ECI to ECEF reference frame for a given time
eci2ecef(double t, const double* rEci, double* rEcef, const double* vEci, double* vEcef);

// Converts positions and velocity from ECEF to ECI reference frame for given time points
// (ld... is the stride; ld = leading dimension, as in FORTRAN)
ecef2eciN(int n, const double* ts, const double* rEcef, int ldEcef, double* rEci, int ldrEci,
	const double* vEcef, int ldvEcef, double* vEci, int ldvEci);

// Converts positions and velocities from ECI to ECEF reference frame for given time points
// (ld... is the stride; ld = leading dimension, as in FORTRAN)
eci2ecefN(int n, const double* ts, const double* rEci, int ldrEci, double* rEcef, int ldEcef,
	const double* vEci, int ldvEci, double* vEcef, int ldvEcef);

// Finds the angular frequency of a satellite with a given semi-major orbital axis
double calcAngFreq(double a);

// Estimates the semi-major orbital axis for a satellite with given trajectory points
double calcMajSemiaxisFromTraj(int n, const double* r3dKm);

// Estimates the angular frequency of a satellite with given trajectory points
double calcAngFreqFromTraj(int n, const double* r3dKm);

// Resamples satellite positions (in ECI) to a new time grid.
// Optionally, finds the velocities (in ECI) at the nodes of the new grid
int resampleSatTrajEci(int nEpochs, const double* ts, double* rEci, int ldrEci,
	int nEpochsOut, const double* tsOut, double* rEciOut, int ldrEciOut, double* vEciOut, int ldvEciOut);

// Resamples satellite positions (in ECEF) to a new time grid.
// Optionally, finds the velocities (in ECEF) at the nodes of the new grid
int resampleSatTrajEcef(int nEpochs, const double* ts, double* rEcef, int ldrEcef,
	int nEpochsOut, const double* tsOut, double* rEcefOut, int ldrEcefOut, double* vEcefOut, int ldvEcefOut);

// Resamples clock values to a new time grid.
// Optionally, finds the clock rates of change at the nodes of the new grid
int resampleSatClock(int nEpochs, const double* ts, double* clock, int ldclock,
	int nEpochsOut, const double* tsOut, double* clockOut, int ldclockOut, double* clockRateOut, int ldclockRateOut);