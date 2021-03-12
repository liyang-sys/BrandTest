#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

int thd_index_inv(int idx)
{
	uint thd_arry[41] = { 0,1,	2,	3,	4,	5	,6	,7,	8,	10,	12,	14,	16,	20,	24,	28,	32,	40,	48,	56,	64,	80,	96,	112,	128,	160,	192	,224	,256	,320	,384	,448	,512	,640,	768	,896,	1024	,1280,	1536,	1792,	2048 };
	return thd_arry[idx];
}
