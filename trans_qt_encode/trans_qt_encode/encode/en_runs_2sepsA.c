#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char


/** @brief
  * @param   r[IN/OUT]
  * @param   p[IN]
  * @param   cbook[IN]
  * @param   sep[OUT]:如果不需要，直接传NULL
  * @retval  none
  * @note    r中数据会被修改，若后面需要，请事先备份
  */
void en_runs_2sepsA(Uint32_Dat* r, float p, unsigned char cbook, SEP* sep)
{
	union data {
		uint a;
		uchar b[4];
	} rem;
	//rem = rem;    // 防止编译warning
	extern int ptr;
	extern uchar *bin;
	int x = 0;
	int thd = 0;
	int cbk = 0;
	int *n = NULL;
	SEP temp;
	GOL gol;
	unsigned int maxr = 0;
	int lenrs = 0;
	int ntem = 0;
	int sum1 = 0, sum2 = 0;
	if (sep == NULL) {
		sep = &temp;
	}

	/*! %separating r into rk and rw ------------------ */
	/* maxr=max(r); */
	maxr = 0;
	for (int i = 0; i <= r->len; i++)
		if (maxr < r->dat[i])   maxr = r->dat[i];
	/* n=hist(r,[1:maxr]); */
	n = hist(r->dat, r->len);
	/* thd=find_thd(n,p); */
	thd = find_thd(r->len, n, p);
	//printf("maxr  = %d thd = %d\n", maxr, thd);
	/*  ntem=max([2*thd-maxr-1, 1]); */
	ntem = 1;
	if ((2 * thd - (int)maxr - 1) > ntem) {
		ntem = 2 * thd - maxr - 1;
	}
	//printf("ntem  = %d \n", ntem);
	/* sum(n(thd:maxr)) */
	sum1 = 0;
	for (int i = thd; i <= maxr; i++) sum1 += n[i];
	/* sum(n(ntem:thd-1)) */
	sum2 = 0;
	for (int i = ntem; i <= thd - 1; i++) sum2 += n[i];
	/* if thd>2 && n(thd)>=n(thd-1) && sum(n(thd:maxr))>=sum(n(ntem:thd-1)) */
	if ((thd > 2) && (n[thd] >= n[thd - 1]) && (sum1 >= sum2))
		/* thd=thd-1; */
		thd += -1;

	/* bin = encode_stationary_source(thd,4);  biny=[biny bin]; */
	encode_stationary_source_bin((unsigned int *)&thd, 1, 4, 0, 0, 0, 0, 0);
	/* [sep,rw,rk]=separate(r,thd);   rw=rw-thd; */
	*sep = separate0((int *)r->dat, thd, r->len, 1);
	/* nk=n(1:thd); */
	Int32_Dat nk;
	nk.len = thd;
	nk.dat = (int *)malloc(sizeof(int)*nk.len);
	memcpy(nk.dat, n + 1, sizeof(int)*nk.len);
	/*! %encoding sep ---------------------------------- */
	/* [rs,lastBit]=Golomb0(sep); */
	gol = Golomb0(sep->sep, sep->lensepbit);
	/* biny=[biny lastBit]; */
	x = ptr & 7;
	bin[ptr >> 3] |= gol.lastBit << (7 - x);
	ptr++;
	/* lenrs=length(rs); */
	lenrs = gol.lenr;
	//printf("lastBit = %d lenrs = %d \n", gol.lastBit, lenrs);
	/* bin = encode_stationary_source(lenrs,cbook);  biny=[biny bin];  */
	/* Hint:此处要求cbook为整数 */
	encode_stationary_source_bin((unsigned int *)&lenrs, 1, cbook, 0, 0, 0, 0, 0);
	/* cbk=4 */
	cbk = 4;
	/* biny=en_runs_1sep(rs,p,cbk,length(sep),biny); */
	Uint32_Dat rs;
	rs.dat = gol.r;
	rs.len = gol.lenr;
	//printf("rs len = %d\n", rs.len);
	en_runs_1sep(&rs, p, cbk, sep->lensepbit, NULL);

	/*! %encoding rk and rw ---------------------------- */
	/* biny=en_runs_sep_sub(thd,rk,nk,rw,biny); */
	en_runs_sep_sub(thd, r->dat, (unsigned int *)nk.dat, sep->rw, sep->lrk, sep->lrw);
#if 0
	内存释放：
		1 - n
		2 - nk.dat
		3 - rs.dat
#endif
		free(n);
	free(nk.dat);
	free(rs.dat);
	if (sep == &temp) {
		free((*sep).rw);
		free((*sep).sep);
	}
	return;
}
