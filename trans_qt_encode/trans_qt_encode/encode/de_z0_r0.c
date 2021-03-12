#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

DEC de_z0_r0(uchar *bin, int len)
{
	extern int ptr;
	DEC decZ;
	DEC decBlock3, decBlock5,dec1;
	decZ.lenzbit = len;
	decZ.z = calloc(((len>>3))+1, sizeof(unsigned char));
	DES des = deSFcode(bin, 36);
	int binBytePos, zBytePos, tempBit, lastBit;
	GOLINV golinv1;
	if (des.sym == 1)
	{
		for (int i = 0; i < len; i++)         //把c移到bin
		{
			binBytePos = ptr & 7;
			zBytePos = i & 7;
			tempBit = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;   //取bin中的一位
			decZ.z[i >> 3] |= (tempBit << (7 - zBytePos));
			ptr++;
		}
	}
	else if(des.sym == 2)
	{
		decBlock3 = decode_block3(bin, len);
		decZ.z = decBlock3.z;
		decZ.lenzbit = decBlock3.lenzbit;
	}
	else if(des.sym == 3)
	{
		decBlock5 = decode_block5(bin, len);
		decZ.z = decBlock5.z;
		decZ.lenzbit = decBlock5.lenzbit;
	}
	else
	{
		binBytePos = ptr & 7;
		lastBit = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;

		dec1 = decode_stationary_source_sumr(bin, des.sym - 4, len);

		golinv1 = GolombInv(dec1.r, lastBit, dec1.lenr);
		decZ.z = golinv1.z;
		decZ.lenzbit = golinv1.lenzbit;
	}
	return decZ;
}