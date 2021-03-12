#include "all.h"
#include <stdio.h>
#include "parameter_setting.h"

#define uint unsigned int
#define uchar unsigned char

Int32_Dat hist_return_len(uint* r, int lenr)
{
	int max = 0;
	for (int i = 0; i < lenr; i++)
	{
		if (max < r[i])
		{
			max = r[i];
		}
	}
	int* n = calloc(max + 1, sizeof(int));
	for (int i = 0; i < lenr; i++)
	{
		n[r[i]]++;
	}
	Int32_Dat ret;
	ret.dat = n;
	ret.len = max + 1;
	return ret;
}

