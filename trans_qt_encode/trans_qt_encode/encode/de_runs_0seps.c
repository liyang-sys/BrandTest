#include <stdio.h>
#include "all.h"
#include "parameter_setting.h"

#define uint unsigned int
#define uchar unsigned char

/** @brief
  * @param   cbook[IN]
  * @param   sep[OUT]
  * @retval  none
  * @note    sep->len存储的是sep->dat的二进制数据长度，
  *          不是字节长度
  */
void de_runs_0seps(unsigned char cbook, Uint8_Dat* sep)
{
	extern uchar *bin;
	extern int  ptr;
	int lastBit = 0;
	int codebook = 0;
	DEC dec;
	DES des;
	GOLINV golinv;
	int lenrs = 0;
	int x1 = 0;

	/* lastBit=biny(ptr); */
	x1 = ptr & 7;
	lastBit = (bin[ptr >> 3] >> (7 - x1)) & 1;
	/* ptr=ptr+1; */
	ptr++;
	/* [lenrs,ptr]=decode_stationary_source_lenr(biny,cbook,1,ptr); */
	dec = decode_stationary_source_lenr(bin, cbook, 1);

	lenrs = dec.r[0];
	free(dec.r);
	dec.r = NULL;
	/* [sym,ptr] = deSFcode(biny,10,ptr); */
	des = deSFcode(bin, 10);
	/* codebook=sym-1; */
	codebook = des.sym - 1;
	/* [rs,ptr]=decode_stationary_source_lenr(biny,codebook,lenrs,ptr); */
	dec = decode_stationary_source_lenr(bin, codebook, lenrs);
	golinv = GolombInv(dec.r, lastBit, dec.lenr);

	// 释放内存
	free(dec.r);

	// 保存结果
	sep->dat = golinv.z;
	sep->len = golinv.lenzbit;

	return;
}
