#include "sp3.h"
#include "error.h"

#pragma warning(disable: 4996)

#define MAX_STR_SIZE 81

static int parseNumEpochs(const char* str)
{
	int nEpoch = 0;
	int _ = sscanf_s(str + 32, "%d", &nEpoch);
	return nEpoch;
}

static int parseNumSats(const char* str)
{
	int nSats = 0;
	int _ = sscanf(str + 3, "%d", &nSats);
	return nSats;
}

static parseTime(const char* str, double* t)
{
	int _ = 0;
	int h = 0;
	int m = 0;
	double s = 0;
	_ = sscanf(str + 14, "%d", &h);
	_ = sscanf(str + 17, "%d", &m);
	_ = sscanf(str + 20, "%lf", &s);
	*t = (h * 60 + m) * 60 + s;
}

static parseData(const char* str, double* data)
{
	double x = 0;
	double y = 0;
	double z = 0;
	double c = 0;
	int _ = 0;
	_ = sscanf(str + 4, "%lf", &x);
	_ = sscanf(str + 18, "%lf", &y);
	_ = sscanf(str + 32, "%lf", &z);
	_ = sscanf(str + 46, "%lf", &c);
	data[0] = x;
	data[1] = y;
	data[2] = z;
	data[3] = c;
}

int loadFromSp3(const char* filename, int* nEpochs, int* nSats, double** ts, double** data)
{
	FILE* fp = fopen(filename, "r");
	CHK_FOPEN(fp, filename);

	char str[MAX_STR_SIZE];
	int lineNum = 0;
	int iEpoch = -1;
	int iSat = -1;
	while (fgets(str, MAX_STR_SIZE, fp) != NULL)
	{
		lineNum++;
		if (lineNum == 1)
		{
			*nEpochs = parseNumEpochs(str);
			*ts = malloc(sizeof(double) * *nEpochs);
			CHK_MALLOC(ts);
		}
		if (lineNum == 3)
		{
			*nSats = parseNumSats(str);
			*data = malloc(sizeof(double) * 4 * *nEpochs * *nSats);
			CHK_MALLOC(data);
		}
		if (str[0] == '*')
		{
			iEpoch++;
			iSat = -1;
			parseTime(str, *ts + iEpoch);
		}
		if (str[0] == 'P')
		{
			iSat++;
			parseData(str, *data + 4 * (iEpoch + *nEpochs * iSat));
		}
	}
	fclose(fp);
	return 0;
}