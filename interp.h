#pragma once

// Finds the value of a piecewise 'lerp' (commonly used name for linear interpolant)
double lerp(int n, const double* x, const double* y, int ldy, double xq);

// Finds the derivative of a piecewise 'lerp' (commonly used name for linear interpolant)
double lerpDer(int n, const double* x, const double* y, int ldy, double xq);

// Finds the value of a trigonometric polynomial
double calcTrigPoly(int m, double omega, const double* a, double t);

// Finds the derivative of a trigonometric polynomial
double calcTrigPolyDer(int m, double omega, const double* a, double t);

// Finds the value of a sliding trigonometric interpolant
double calcSlidingTrigInterp(int n, int k, double omega, const double* a, const double* ts, double t);

// Finds the derivative of a sliding trigonometric interpolant
double calcSlidingTrigInterpDer(int n, int k, double omega, const double* a, const double* ts, double t);

// Finds the value of a sliding trigonometric interpolant
int findTrigPoly(int m, double omega, const double* t, const double* p, int ldp, double* a);

// Finds the derivative of a sliding trigonometric interpolant
int findSlidingTrigInterp(int n, int k, double omega, const double* t, const double* p, int ldp, double* a);