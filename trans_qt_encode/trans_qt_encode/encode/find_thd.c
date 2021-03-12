#include "all.h"
#include <stdio.h>

#define uint unsigned int
#define uchar unsigned char

int find_thd (int lenr, int* n, float p)//
{
	float psm = p * lenr;
	int thd = 1;
	int ssm = n[1];
	while (ssm < psm)
	{
		thd++;
		ssm += n[thd];
	}
	return thd;
}

