#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DEC de_r1 (uchar *biny, int lenbinbit)
{
    extern int ptr;
    int x = 0; ///一个字节里当前的位数
    int i;
    int lenrk;
    int lenrw;
    int thd;
    int idx;
    int lastBit;
    int codebook;
    uint *r;
    DES des1;
    DEC dec1, dec2;
    GOLINV golinv1;
    DRSS drss1;

    des1 = deSFcode (biny, 40);
    idx = des1.sym;
    thd = thd_index_inv (idx);

    /************ decoding sep ************/

    x = ptr & 7;
    lastBit = ((biny[ptr >> 3] >> (7 - x)) & 1);
    ptr++;
    dec1 = decode_stationary_source_lenr (biny, 20, 1); ///dec1.r[0] = lenrs
    des1 = deSFcode (biny, 5);
    des1.sym--; ///des1.sym = codebook
    dec2 = decode_stationary_source_lenr (biny, des1.sym, dec1.r[0]); ///dec2.r = rs, dec2.lenr = lenrs
    golinv1 = GolombInv (dec2.r, lastBit, dec2.lenr); ///golinv1.z = sep, golinv1.lenzbit = lensepbit

    /************ rw and rk ************/

    lenrw = dec2.lenr - (!lastBit); ///dec2.lenr = lenrs
    lenrk = golinv1.lenzbit - lenrw; ///golinv1.lenzbit = lensepbit
    drss1 = de_runs_sep_sub (biny, thd, lenrk, lenrw, lenbinbit);

    /************ synthesis r ************/

    r = separate_inv (golinv1.z, drss1.rw, drss1.rk, golinv1.lenzbit, thd); ///golinv1.z = sep, golinv1.lenzbit = lensepbit

	//添加的
	DEC ret;
	ret.r = r;
	ret.lenr = golinv1.lenzbit;
    return ret;
}
