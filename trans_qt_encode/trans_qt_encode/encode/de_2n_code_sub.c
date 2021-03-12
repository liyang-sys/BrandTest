#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
#include <math.h>

DEC de_2n_code_sub(uchar *bin, int n,Uint8_Dat* z)
{
	extern int ptr;
	int lenw = 0;
	int zBytePos = 0;
	int binBytePos = ptr & 7;
	for (int i = 0; i < z->len; i++)
	{
		zBytePos = i & 7;
		int temp = ((z->dat[i >> 3]) >> (7 - zBytePos)) & 1;
		if (temp == 1)
		{
			lenw++;
		}
	}
	int lenk = z->len - lenw;

	Uint32_Dat cw;
	cw.len = lenw;
	cw.dat = calloc(lenw, sizeof(unsigned int));

	Uint32_Dat ck;
	ck.len = lenk;
	ck.dat = calloc(lenk, sizeof(unsigned int));
	unsigned int temp = 0;
	for (int i = 0; i < lenw; i++)
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
		cw.dat[i] = temp + 1;
	}
	for (int i = 0; i < lenk; i++)
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
		ck.dat[i] = temp +1;
	}
	
	DEC c;
	c.lenr = z->len;
	c.r = separate_inv(z->dat, cw.dat, ck.dat, c.lenr,(int)pow(2,n));
	if (cw.dat != NULL)
	{
		free(cw.dat);
	}
	if (ck.dat != NULL)
	{
		free(ck.dat);
	}
	return c;

}