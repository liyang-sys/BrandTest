#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char


/** @brief
  * @param   r[IN/OUT]
  * @param   p[IN]
  * @param   cbook[IN]
  * @param   lenz[IN]
  * @param   sep[OUT]:如果不需要，直接传NULL
  * @retval  none
  * @note    r中数据会被修改，若后面需要，请事先备份
  */
void en_runs_2seps(Uint32_Dat* r, float p, int lenz, unsigned char cbook, SEP* sep)
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
	int lenrs = 0;
	if (sep == NULL) {
		sep = &temp;
	}

	/*! %if lenz > 1000 */
	if (lenz > 2000) {
		/* biny=[biny 0]; */
		ptr++;
		/*! %separating r into rk and rw ------------------ */
		/* n=hist(r,[1:max(r)]); */
		n = hist(r->dat, r->len);
		/* thd=find_thd(n,p); */
		thd = find_thd(r->len, n, p);
		/* bin = encode_stationary_source(thd,4);  biny=[biny bin]; */
		encode_stationary_source_bin((unsigned int *)&thd, 1, 4, 0, 0, 0, 0, 0);
		/* [sep,rw,rk]=separate(r,thd); + rw=rw-thd; */
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
		/* bin = encode_stationary_source(lenrs,cbook);  biny=[biny bin]; */
		/* Hint:此处要求cbook为整数 */
		encode_stationary_source_bin((unsigned int *)&lenrs, 1, cbook, 0, 0, 0, 0, 0);
		/* cbk = 4 */
		cbk = 4;
		/* biny=en_runs_1sep(rs,p,cbk,length(sep),biny); */
		Uint32_Dat rs;
		rs.dat = gol.r;
		rs.len = gol.lenr;
		en_runs_1sep(&rs, p, cbk, sep->lensepbit, NULL);

		/*! biny=en_runs_sep_sub(thd,rk,nk,rw,biny);  */
		en_runs_sep_sub(thd, r->dat, (unsigned int *)nk.dat, sep->rw, sep->lrk, sep->lrw);
#if 0
		内存释放:
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
	}
	else {
		/* biny=[biny 1]; */
		x = ptr & 7;
		bin[ptr >> 3] |= 1 << (7 - x);
		ptr++;
		/* cbk=4; */
		cbk = 4;
		/* biny=en_runs_1sep(r,p,cbk,lenz,biny); */
		en_runs_1sep(r, p, cbk, lenz, NULL);
	}


	return;
}
