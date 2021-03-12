#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
#include <math.h>

DEC de_2n_code_sub1(uchar *bin, int n, int len)
{
	extern int ptr;
	int binBytePos = ptr & 7;
	DEC cr;
	cr.lenr = len;
	cr.r = calloc(len, sizeof(unsigned int));
	int temp = 0;
	for (int i = 0; i < len; i++)
	{
		temp = 0;
		for (int j = 0; j < n; j++)
		{
			binBytePos &= 7;
			temp = temp << 1;
			temp += (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
			ptr++;
			binBytePos++;
		}
		cr.r[i] = temp + 1;
	}
	return cr;
}