#include <stdio.h>
#include "all.h"
#include "parameter_setting.h"
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
void de_runs_1sep(int lenr, unsigned char cbook, Uint32_Dat* r, int lenbinbit)
{
	extern uchar *bin;
	extern int  ptr;
	r->len = lenr;
	Uint32_Dat  thd;
	Uint32_Dat  rs;
	DEC dec;
	DES des;
	DRSS drss;
	int lenrw = 0, lenrk = 0;
	GOLINV sep;
	int lastBit = 0;
	unsigned int codebook = 0;

	int flg = 0;
	unsigned char x1;
	/* flg=biny(ptr); */
	x1 = ptr & 7;
	flg = (bin[ptr >> 3] >> (7 - x1)) & 1;
	/* ptr=ptr+1; */
	ptr++;


	if (flg < 0.5) {   /* %length(z)>600 */
	   /* [thd,ptr]=decode_stationary_source_lenr(biny,cbook,1,ptr); */
		dec = decode_stationary_source_lenr(bin, cbook, 1);
		thd.len = dec.lenr;
		thd.dat = dec.r;
		dec.r = NULL;
		/* %decoding sep -------------------------------------------------------- */
		/* lastBit=biny(ptr); */
		x1 = ptr & 7;
		lastBit = (bin[ptr >> 3] >> (7 - x1)) & 1;
		/* ptr=ptr+1; */
		ptr++;
		/* [codebook,ptr] = deSFcode(biny,5,ptr); */
		des = deSFcode(bin, 5);
		codebook = des.sym;
		/* codebook=codebook-1; */
		codebook += -1;

		/* [rs,ptr]=decode_stationary_source_sumr(biny,codebook,lenr,ptr); */
		dec = decode_stationary_source_sumr(bin, codebook, lenr);
		rs.len = dec.lenr;
		rs.dat = dec.r;
		dec.r = NULL;
		/* sep=GolombInv(rs,lastBit); */
		sep = GolombInv(rs.dat, lastBit, rs.len);

		/* %decoding rw and rk -------------------------------------------------- */
		/* [rk,rw,ptr]=de_runs_sep_sub(biny,ptr,thd,sep); */
		lenrw = rs.len - (!lastBit);
		lenrk = sep.lenzbit - lenrw;
		drss = de_runs_sep_sub(bin, thd.dat[0], lenrk, lenrw, lenbinbit);

		/*! %synthesis r -------------------------------------------------------- */
		/* r=separate_inv(sep,rw+thd,rk); */
		r->dat = separate_inv(sep.z, (int *)drss.rw, (int *)drss.rk, lenr, thd.dat[0]);

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
	else {   /* %length(z)<600 */
	  /* [codebook,ptr] = deSFcode(biny,6,ptr); */
		des = deSFcode(bin, 6);
		codebook = des.sym;
		/* codebook=codebook-1; */
		codebook += -1;
		/* [r,ptr]=decode_stationary_source_lenr(biny,codebook,lenr,ptr); */
		dec = decode_stationary_source_lenr(bin, codebook, lenr);
		r->len = dec.lenr;
		r->dat = dec.r;
		dec.r = NULL;
	}

	return;
}
