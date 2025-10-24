#include "sp3.h"
#include "astrodyn.h"
#include "error.h"
#include <stdio.h>

#pragma warning(disable: 4996)

static linspace(double a, double b, int n, double* x)
{
	if (n <= 0)
	{
		return 0;
	}
	if (n == 1)
	{
		x[0] = a;
		return 0;
	}
	double h = (b - a) / (n - 1);
	x[0] = a;
	x[n - 1] = b;
	for (int i = 1; i < n - 1; i++)
	{
		x[i] = a + i * h;
	}
	return 0;
}

static saveNumbers(const char* filename, int n, const double* a)
{
	FILE* fp = fopen(filename, "w");
	CHK_FOPEN(fp, filename);

	for (int i = 0; i < n; i++)
	{
		fprintf(fp, "%13lf\n", a[i]);
	}

	fclose(fp);
	return 0;
}

static saveVectors3d(const char* filename, int n, const double* v)
{
	FILE* fp = fopen(filename, "w");
	CHK_FOPEN(fp, filename);

	for (int i = 0; i < n; i++)
	{
		const double* vi = v + 3 * i;
		double vx = vi[0];
		double vy = vi[1];
		double vz = vi[2];
		fprintf(fp, "%13lf %13lf %13lf\n", vx, vy, vz);
	}

	fclose(fp);
	return 0;
}

int main(int argc, char* argv[])
{
	const char srcFilenameDefault[] = "../../../test.sp3";
	const int satDefault = 2;
	const char posFilenameDefault[] = "positions.txt";
	const char velFilenameDefault[] = "velocities.txt";
	const char clockFilenameDefault[] = "clock.txt";
	const char clockRateFilenameDefault[] = "clockRate.txt";

	const char* srcFilename = argc > 1 ? argv[1] : srcFilenameDefault;
	int sat = argc > 2 ? strtol(argv[2], NULL, 10) : satDefault;
	const char* posFilename = argc > 3 ? argv[3] : posFilenameDefault;
	const char* velFilename = argc > 4 ? argv[4] : velFilenameDefault;
	const char* clockFilename = argc > 5 ? argv[5] : clockFilenameDefault;
	const char* clockRateFilename = argc > 6 ? argv[6] : clockRateFilenameDefault;

	int nEpochs = 0;
	int nSats = 0;
	double* ts = NULL;
	double* data = NULL;
	if (loadFromSp3(srcFilename, &nEpochs, &nSats, &ts, &data))
	{
		return 1;
	}
	printf("Loaded data for %d epochs\nof satellite #%d\nfrom '%s'\n",
		nEpochs, sat, srcFilename);

	double* r = data + 4 * nEpochs * sat;
	double* clock = r + 3;

	int nEpochsOut = 3 * nEpochs - 2;
	double* tsOut = (double*)malloc(sizeof(double) * nEpochsOut);
	CHK_MALLOC(tsOut);

	linspace(0, ts[nEpochs - 1], nEpochsOut, tsOut);

	double* rOut = (double*)malloc(sizeof(double) * 3 * nEpochsOut);
	CHK_MALLOC(rOut);

	double* vOut = (double*)malloc(sizeof(double) * 3 * nEpochsOut);
	CHK_MALLOC(vOut);

	double* clockOut = (double*)malloc(sizeof(double) * nEpochsOut);
	CHK_MALLOC(clockOut);

	double* clockRateOut = (double*)malloc(sizeof(double) * nEpochsOut);
	CHK_MALLOC(clockRateOut);

	printf("\n");

	resampleSatTrajEcef(nEpochs, ts, r, 4, nEpochsOut, tsOut, rOut, 3, vOut, 3);
	printf("Resampled satellite positions.\nFound velocities\n");

	printf("\n");

	resampleSatClock(nEpochs, ts, clock, 4, nEpochsOut, tsOut, clockOut, 1, clockRateOut, 1);
	printf("Resampled clock values.\nFound clock rates of change\n");

	printf("\n");

	saveVectors3d(posFilename, nEpochsOut, rOut);
	printf("Saved positions to '%s'\n", posFilename);
	
	saveVectors3d(velFilename, nEpochsOut, vOut);
	printf("Saved velocities to '%s'\n", velFilename);

	saveNumbers(clockFilename, nEpochsOut, clockOut);
	printf("Saved clock values to '%s'\n", clockFilename);

	saveNumbers(clockRateFilename, nEpochsOut, clockRateOut);
	printf("Saved clock rates of change to '%s'\n", clockRateFilename);

	free(ts);
	free(data);
	free(rOut);
	free(vOut);
	free(clockOut);
	free(clockRateOut);
	return 0;
}