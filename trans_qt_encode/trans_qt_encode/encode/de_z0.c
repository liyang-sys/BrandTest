#include <stdio.h>
#include "all.h"
#include "encoding.h"

DEC de_z0(uchar *bin, int len)
{
	extern int ptr;
	//printf("\n de_z0的 ptr = %d\n", ptr);
	DEC decZ;
	DEC decBlock3, decBlock5,dec1;
	GOLINV golinv1;
	decZ.z = calloc((len), sizeof(int));
	int binBytePos = ptr & 7;
	int zBytePos = 0;
	int tempBit,bn,lastBit;

	bn = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;    //取bin中的某一位(cbook的高位)
	bn = bn << 1;
	ptr++;

	binBytePos = ptr & 7;
	bn = ((bin[ptr >> 3] >> (7 - binBytePos)) & 1) + bn;  //(cbook第二位，加一
	bn = bn << 1;
	ptr++;

	binBytePos = ptr & 7;
	bn = ((bin[ptr >> 3] >> (7 - binBytePos)) & 1) + bn;  //(cbook的低位，加一
	ptr++;

	switch (bn)
	{
	case 0:
		for (int i = 0; i < len; i++)         //把c移到bin
		{
			binBytePos = ptr & 7;
			zBytePos = i & 7;
			tempBit = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;   //取bin中的一位
			decZ.z[i >> 3] |= (tempBit << (7 - zBytePos));
			ptr++;
		}
		decZ.lenzbit = len;
		break;
	case 1:
		decBlock3 = decode_block3(bin, len);
		decZ.z = decBlock3.z;
		decZ.lenzbit = decBlock3.lenzbit;
		break;
	case 2:
		decBlock5 = decode_block5(bin, len);
		decZ.z = decBlock5.z;
		decZ.lenzbit = decBlock5.lenzbit;
		break;
	case 3:
		binBytePos = ptr & 7;
		lastBit = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;

		dec1 = decode_stationary_source_sumr(bin, 0, len);

		golinv1 = GolombInv(dec1.r, lastBit, dec1.lenr);
		decZ.z = golinv1.z;
		decZ.lenzbit = golinv1.lenzbit;
		break;
	case 4:
		binBytePos = ptr & 7;
		lastBit = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;

		dec1 = decode_stationary_source_sumr(bin, 1, len);

		golinv1 = GolombInv(dec1.r, lastBit, dec1.lenr);
		decZ.z = golinv1.z;
		decZ.lenzbit = golinv1.lenzbit;
		break;
	case 5:
		binBytePos = ptr & 7;
		lastBit = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;

		dec1 = decode_stationary_source_sumr(bin, 2, len);

		golinv1 = GolombInv(dec1.r, lastBit, dec1.lenr);
		decZ.z = golinv1.z;
		decZ.lenzbit = golinv1.lenzbit;
		break;
	case 6:
		binBytePos = ptr & 7;
		lastBit = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;

		dec1 = decode_stationary_source_sumr(bin, 3, len);

		golinv1 = GolombInv(dec1.r, lastBit, dec1.lenr);
		decZ.z = golinv1.z;
		decZ.lenzbit = golinv1.lenzbit;
		break;
	default:
		binBytePos = ptr & 7;
		lastBit = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;

		dec1 = decode_stationary_source_sumr(bin, 4, len);

		golinv1 = GolombInv(dec1.r, lastBit, dec1.lenr);
		decZ.z = golinv1.z;
		decZ.lenzbit = golinv1.lenzbit;
		break;
	}
	//printf("\n 解码之后的 ptr = %d\n", ptr);

	binBytePos = ptr & 7;
	int inv = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
	ptr++;
	if (inv)
	{
		//printf("\ninv长度 %d\n", decZ.lenzbit >> 3);
		for (int i = 0; i <= decZ.lenzbit >> 3; i++)   
		{
			decZ.z[i] = ~decZ.z[i];
		}
	}

	return decZ;
}