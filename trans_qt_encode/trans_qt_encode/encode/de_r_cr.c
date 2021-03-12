#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

DEC de_r_cr(uchar *bin, int len,int lenbinbit)
{
	DEC ret;
	ret.z = NULL;
	ret.lenzbit = 0;
	extern int ptr;
	int binBytePos = ptr & 7;
	int flg = 0;
	flg = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;  
	flg = flg << 1;
	ptr++;

	binBytePos = ptr & 7;
	flg = ((bin[ptr >> 3] >> (7 - binBytePos)) & 1) + flg; 
	ptr++;
	if (flg == 3)
	{
		ret = de_zone_new(bin,lenbinbit);
	}
	else if (flg == 2)
	{
		binBytePos = ptr & 7;
		int lB = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;
		ret = de_zone_sub2(bin, lB, lenbinbit);
	}
	else if (flg == 1)
	{
		binBytePos = ptr & 7;
		int lB = 0;
		lB = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;

		DEC dec = de_r1(bin, lenbinbit);
		GOLINV gol = GolombInv(dec.r, lB, dec.lenr);
		ret.z = gol.z;
		ret.lenzbit = gol.lenzbit;
		if (dec.r != NULL && dec.lenr != 0)
		{
			free(dec.r);
		}
	}
	else
	{
		ret = de_z0_r0(bin, len);
	}
	return ret;
}