#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

/** @brief
  * @param   r
  * @param   p :
  * @param   lenr : r数组的长度
  * @param   n:      r的频率分布直方图
  * @param   on_off：1时不做hist，0时做hist
  * @return
  * @note
  */
int en_r1 (uint *r, float p, int lenr, int *n, int on_off)
{
    union data {
      unsigned short int a;
      uchar b[4];
    } rem;
    extern int ptr;
    extern uchar *bin;
    int i;
    int x; ///编码时：一个字节里当前的位数
    int thd;
    int idx;
    GOL gol1;
    SFC sfc1;
    THE the1;
    SEP sep1;
    ESSC essc1;

    /************ separating r into rk and rw ************/

    if (!on_off)
    {
        n = hist (r, lenr);
    }
	thd = find_thd (lenr, n, p);
	the1 = the_index (thd);
	thd = the1.thd;
	idx = the1.idx;
	sfc1 = SFcode (idx, 40);
    x = ptr & 7;


    /* 将这段编码移动到最前面 */
    rem.a = sfc1.code;
    rem.a = rem.a << (16 - x - sfc1.lb);

    bin[ptr >> 3] |= rem.b[1];
    bin[(ptr >> 3) + 1] |= rem.b[0];

    /* 指针偏移，位数偏移，位数限制 */
    ptr += sfc1.lb;
    x += sfc1.lb;
    x &= 7;

	sep1 = separate0 (r, thd, lenr, 1); ///rk存放在r里面


	uint *nk = calloc (thd, sizeof(uint));
	for (i = 0; i < thd; i++)
	{
		nk[i] = n[i + 1];
	}

    /************ encoding sep ************/

    gol1 = Golomb0 (sep1.sep, sep1.lensepbit); ///lenrs是gol1.lenr，rs是 gol1.r, ps是gol1.p
    bin[ptr >> 3] |= gol1.lastBit << (7 - x);
    ptr++; x++; x &= 7;

    encode_stationary_source_bin (&gol1.lenr, 1, 20, 0, 0, 0, 0, 0); ///biny=[biny bin_lenrs]
    essc1 = encode_stationary_source_cbook (gol1.p);
    sfc1 = SFcode ((essc1.codebook + 1), 5);
    x = ptr & 7;
    rem.a = sfc1.code; ///biny=[biny codebook]
    rem.a = rem.a << (16 - x - sfc1.lb);
    bin[ptr >> 3] |= rem.b[1];
    bin[(ptr >> 3) + 1] |= rem.b[0];
    ptr += sfc1.lb; x += sfc1.lb; x &= 7;
    encode_stationary_source_bin (gol1.r, gol1.lenr, gol1.p, essc1.k, essc1.m, essc1.m1, essc1.m2, essc1.cls); ///biny=[biny rs]

    /************ encoding rk and rw ************/

    en_runs_sep_sub (thd, r, nk, sep1.rw, sep1.lrk, sep1.lrw);
}

