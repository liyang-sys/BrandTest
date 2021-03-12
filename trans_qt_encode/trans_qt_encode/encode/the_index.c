#include "all.h"
#include <stdio.h>

#define uint unsigned int
#define uchar unsigned char

THE the_index(int thd)
{
    int idx;
	THE the1 = { 0 };
	int thd_arry[41] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 320, 384, 448, 512, 640, 768, 896, 1024, 1280, 1536, 1792, 2048};


	if (thd > 2048)
	{
		thd = 2048;
	}
	idx = 1;
	while (thd_arry[idx] < thd)
	{
		idx++;
	}
	if (idx > 6)
	{
		if ((thd_arry[idx] - thd) > (thd - thd_arry[idx - 1]))
		{
			idx--;
		}
	}
	the1.thd = thd_arry[idx];
	the1.idx = idx;
	return the1;
}
