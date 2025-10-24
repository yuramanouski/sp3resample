#pragma once
#include <stdio.h>

#define CHK_MALLOC(p) if ((p) == NULL) { \
	printf("Error in file %s, line %d: could not allocate memory", __FILE__, __LINE__); \
	return 1; }

#define CHK_FOPEN(fp, filename) if ((fp) == NULL) { \
	printf("Error in file %s, line %d: could not open file %s", __FILE__, __LINE__, filename); \
	return 1; }
