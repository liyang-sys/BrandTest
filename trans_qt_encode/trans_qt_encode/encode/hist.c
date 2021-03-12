#include "all.h"
#include <stdio.h>

#define uint unsigned int
#define uchar unsigned char

int* hist(uint* r,int lenr) //位置零没有用，得从1开始
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

	return n;
}

