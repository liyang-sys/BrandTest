#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char

/** @brief
  * @param   lenr[IN]
  * @param   cbook[IN]
  * @param   r[OUT]
  * @param   lenbinbit[IN]:bin的二进制总长度
  * @retval  none
  * @note
  */
void de_runs_2seps(int lenr, unsigned char cbook, Uint32_Dat* r, int lenbinbit)
{
	extern uchar *bin;
	extern int  ptr;
	DEC dec;
	Uint32_Dat  thd;
	DRSS drss;
	DES des;
	int lastBit = 0;
	int lenrs = 0;
	int flg = 0;
	unsigned char x1;
	unsigned char cbk = 0;
	int lenrw = 0, lenrk = 0;
	unsigned char codebook = 0;
	GOLINV sep;
	/* flg=biny(ptr); */
	x1 = ptr & 7;
	flg = (bin[ptr >> 3] >> (7 - x1)) & 1;
	/* ptr=ptr+1; */
	ptr++;
	r->len = lenr;
	if (flg < 0.5) {  /* %length(z)>1000 */
	   /* [thd,ptr]=decode_stationary_source_lenr(biny,4,1,ptr); */
		dec = decode_stationary_source_lenr(bin, 4, 1);
		thd.len = dec.lenr;
		thd.dat = dec.r;
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
		drss = de_runs_sep_sub(bin, thd.dat[0], lenrk, lenrw, lenbinbit);

		/*! %synthesis r -------------------------------------------------------- */
		/* r=separate_inv(sep,rw+thd,rk); */
		r->dat = separate_inv(sep.z, (int *)drss.rw, (int *)drss.rk, lenr, thd.dat[0]);

		/* 内存释放 */
		free(thd.dat);
		free(rs.dat);
		free(sep.z);
		free(drss.rw);
		free(drss.rk);
#if 0
		内存释放：
			1 - thd.dat
			2 - rs.dat
			3 - sep.z
			4 - drss.rw
			5 - drss.rk
#endif
	}
	else {   /* %length(z)<1000 */
	  /* cbk=4; */
		cbk = 4;
		/* [r,ptr]=de_runs_1sep(biny,lenr,cbk,ptr); */
		//de_runs_1sep(lenr, cbook, r, lenbinbit);
		de_runs_1sep(lenr, cbk, r, lenbinbit);
	}

	return;
}
