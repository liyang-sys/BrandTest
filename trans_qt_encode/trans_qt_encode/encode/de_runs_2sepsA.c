#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char

/** @brief
  * @param   cbook[IN]
  * @param   r[OUT]
  * @param   lenbinbit[IN]:bin的二进制总长度
  * @retval  none
  * @note
  */
void de_runs_2sepsA(unsigned char cbook, Uint32_Dat* r, int lenbinbit)
{
	extern uchar* bin;
	extern int  ptr;
	DEC dec;
	int  thd = 0;
	DRSS drss;
	int lastBit = 0;
	int lenrs = 0;
	unsigned char x1;
	unsigned char cbk = 0;
	int lenrw = 0, lenrk = 0;
	GOLINV sep;

	/* [thd,ptr]=decode_stationary_source_lenr(biny,4,1,ptr); */
	dec = decode_stationary_source_lenr(bin, 4, 1);
	thd = dec.r[0];
	free(dec.r);
	dec.r = NULL;

	/*! %decoding sep -------------------------------------------------------- */
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
	/* cbk=4; */
	cbk = 4;
	/* de_runs_1sep(biny,lenrs,cbk,ptr);  */
	Uint32_Dat rs;
	rs.len = lenrs;
	de_runs_1sep(rs.len, cbk, &rs, lenbinbit);
	/* sep=GolombInv(rs,lastBit); */
	sep = GolombInv(rs.dat, lastBit, rs.len);

	/*! %decoding rw and rk -------------------------------------------------- */
	/* [rk,rw,ptr]=de_runs_sep_sub(biny,ptr,thd,sep); */
	lenrw = rs.len - (!lastBit);
	lenrk = sep.lenzbit - lenrw;
	//printf("len rw = %d len rk = %d \n", lenrw, lenrk);
	drss = de_runs_sep_sub(bin, thd, lenrk, lenrw, lenbinbit);

	/*! %synthesis r -------------------------------------------------------- */
	/* r=separate_inv(sep,rw+thd,rk); */
	r->dat = separate_inv(sep.z, (int *)drss.rw, (int *)drss.rk, sep.lenzbit, thd);
	r->len = sep.lenzbit;
	//printf("len r =%d\n", r->len);

	/* 内存释放 */
	free(rs.dat);
	free(sep.z);
	free(drss.rw);
	free(drss.rk);

#if 0
	内存释放：
		1 - rs.dat
		2 - sep.z
		3 - drss.rw
		4 - drss.rk
#endif

		return;
}
