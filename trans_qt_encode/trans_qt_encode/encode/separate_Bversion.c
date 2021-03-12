#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

SEP separate_Bversion(uint *r,uint thd,int r_len, uchar* sep, uint* rw, uint* rk)//在separate中就减去thd
{
	int ptrw = 0, ptrk = 0;
	int  index = 0;
	for (int i = 0; i < r_len; i++)
	{
		if (r[i] > thd)
		{
			sep[i>>3]|=1<<(7-index);
			index++;
			index &= 7;
			rw[ptrw] = r[i];
			rw[ptrw] -= thd;
			ptrw++;
		}
		else
		{
			rk[ptrk] = r[i];
			ptrk++;
			index++;
			index &= 7;
		}
	}
	SEP sep1 = { 0 };
	sep1.lrw = ptrw;
	sep1.lrk= ptrk;
	return sep1;
}
