#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

uint *separate_inv(uchar *sep, int *rw, int *rk, int len, int thd)
{
   int i, t, t0;
   int pw = 0, pk = 0;

   uint *r = (uint *)calloc(len, sizeof(uint));
   //uint *r = (uint *)malloc(len * sizeof(uint));
   int m,n;
   for (i = 0; i < len; i++)
   {
      t0 = i & 7;
      t = (sep[i >> 3] >> (7 - t0)) & 1;
      if (!t){
         r[i] = rk[pk];
         pk++;
      }
      else{
         r[i] = rw[pw] + thd;
         pw++;
      }
   }
   return r;
}

/*
uint *separate_inv (uchar *sep, uint *rw, uint *rk, int lensepbit)
{
	int ptrk = 0, ptrw = 0, i, temp, index = 0;
	uint *r;
	r = (uint*)calloc(lensepbit, sizeof(uint));
	if (!r){
		printf("创建数组r失败！\n");
		exit(1);
	}
	for (i = 0; i < lensepbit; i++){
		temp = (sep[i >> 3] >> (7 - index)) & 0x01;
		if (temp < 0.5){
			r[i] = rk[ptrk];
			ptrk++;
		}
		else{
			r[i] = rw[ptrw];
			ptrw++;
		}
		index++;
		index &= 7;
	}
	return r;
}
*/
