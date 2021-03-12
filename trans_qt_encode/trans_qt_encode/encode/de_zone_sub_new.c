#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char

static unsigned char Table_sep[256] = {
									   0, 1, 1, 2, 1, 2, 2, 3,
									   1, 2, 2, 3, 2, 3, 3, 4,
									   1, 2, 2, 3, 2, 3, 3, 4,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   1, 2, 2, 3, 2, 3, 3, 4,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   1, 2, 2, 3, 2, 3, 3, 4,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   4, 5, 5, 6, 5, 6, 6, 7,
									   1, 2, 2, 3, 2, 3, 3, 4,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   4, 5, 5, 6, 5, 6, 6, 7,
									   2, 3, 3, 4, 3, 4, 4, 5,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   4, 5, 5, 6, 5, 6, 6, 7,
									   3, 4, 4, 5, 4, 5, 5, 6,
									   4, 5, 5, 6, 5, 6, 6, 7,
									   4, 5, 5, 6, 5, 6, 6, 7,
									   5, 6, 6, 7, 6, 7, 7, 8
};

/** @brief
  * @param  _bin[IN]
  * @param  lenbinbit[IN]:bin二进制长度
  * @param  rrw[OUT]
  * @param  rk[OUT]
  * @retval
  * @note
  */
void de_zone_sub_new(uchar *_bin, int lenbinbit, Uint32_Dat *rrw, Uint32_Dat *rk)
{
	extern int ptr;
	extern unsigned char* bin;
	bin = _bin;

	DES  des;
	DEC  dec;
	DER  der;
	GOLINV golinv;
	int nsep = 0;
	int *thd = NULL;
	Uint8_Dat* rsep = NULL;
	unsigned char cbook = 0;
	unsigned char codebook = 0;
	unsigned char lB1 = 0;
	unsigned char lB2 = 0;
	int lenr1 = 0;
	int lenr2 = 0;
	int x = 0;
	int sum = 0;
	Uint32_Dat rs1;
	Uint32_Dat rs2;
	int *lenk = NULL;
	int flg = 0;
	int th[25] = { 0,1,2,3,4,6,8,12,16,24,32,48,64,
				  96,128,192,256,384,512,768,1024,
				  1536,2048,3072,4096 };
	des = deSFcode(_bin, 5);
	nsep = des.sym + 1;
	thd = malloc(sizeof(int)*(nsep + 1));
	if (thd == NULL) {
		perror("thd");
		exit(-1);
	}
	memset(thd, 0, sizeof(int)*(nsep + 1));

	for (int nsp = 1; nsp <= nsep; nsp++) {
		dec = decode_stationary_source_lenr(_bin, 3, 1);
		thd[nsp] = th[dec.r[0]];
	}

	/* rsep=cell(1,nsep); */
	rsep = malloc(sizeof(Uint8_Dat)*nsep);
	memset(rsep, 0, sizeof(Uint8_Dat)*nsep);
	/*! %decode rsep{nsep}: ------------------------------------------------ */
	cbook = 17;

	de_runs_0seps(cbook, &(rsep[nsep - 1]));

	/*! %decode rsep{k} (k from 2 to nsep-1): */
	for (int nsp = nsep - 1; nsp >= 2; nsp--) {
		x = ptr & 7;
		lB2 = (_bin[ptr >> 3] >> (7 - x)) & 1;
		ptr++;
		lenr2 = rsep[nsp].len + (1 - lB2);
		cbook = 4;
		de_runs_1sep(lenr2, cbook, &rs2, lenbinbit);
		golinv = GolombInv(rs2.dat, lB2, rs2.len);
		rsep[nsp - 1].dat = golinv.z;
		rsep[nsp - 1].len = golinv.lenzbit;
		/* 内存释放 */
		free(rs2.dat);
		rs2.dat = NULL;
		rs2.len = 0;
	}
	/*! %decode rsep{1}: --------------------------------------------------- */
	x = ptr & 7;
	lB1 = (_bin[ptr >> 3] >> (7 - x)) & 1;
	ptr++;
	lenr1 = rsep[1].len + (1 - lB1);
	cbook = 17;

	de_runs_2seps(lenr1, cbook, &rs1, lenbinbit);
	golinv = GolombInv(rs1.dat, lB1, rs1.len);
	rsep[0].dat = golinv.z;
	rsep[0].len = golinv.lenzbit;

	/* 内存释放 */
	free(rs1.dat);
	rs1.dat = NULL;
	rs1.len = 0;

	/*! %decode rw: -------------------------------------------------------- */
	/* sum(rsep{nsep}) */
	sum = 0;
	for (int i = 0; i < (rsep[nsep - 1].len + 7) / 8; i++) {
		sum += Table_sep[rsep[nsep - 1].dat[i]];
	}

	der = de_r0(_bin, sum, lenbinbit);
	rrw->dat = der.r;
	rrw->len = der.lenr;

	/* %decoding rk's: ---------------------------------------------------- */
	lenk = malloc(sizeof(int)*nsep);
	memset(lenk, 0, sizeof(int)*nsep);
	for (int nsp = 1; nsp <= nsep; nsp++) {
		sum = 0;
		for (int i = 0; i < (rsep[nsp - 1].len + 7) / 8; i++) {
			sum += Table_sep[rsep[nsp - 1].dat[i]];
		}
		lenk[nsp - 1] = rsep[nsp - 1].len - sum;
	}

	rk = malloc(sizeof(Uint32_Dat)*nsep);
	if (rk == NULL) {
		perror("rk");
		exit(-1);
	}
	/*! %decode rk{3} and rk{2} */
	for (int nsp = nsep; nsp >= 1; nsp--) {
		if (thd[nsp] > 1) {
			x = ptr & 7;
			flg = (_bin[ptr >> 3] >> (7 - x)) & 1;
			ptr++;
			if (flg == 0) {
				dec = de_Kside_new(_bin, lenk[nsp - 1], thd[nsp], lenbinbit);
				rk[nsp - 1].dat = dec.r;
				rk[nsp - 1].len = dec.lenr;
				dec.r = NULL;
				dec.lenr = 0;
			}
			else {
				des = deSFcode(_bin, 20);
				codebook = des.sym;
				codebook += -1;
				dec = decode_stationary_source_Nsym_lenr(_bin, codebook, lenk[nsp - 1], thd[nsp], lenbinbit);
				rk[nsp - 1].dat = dec.r;
				rk[nsp - 1].len = dec.lenr;
				dec.r = NULL;
				dec.lenr = 0;
			}
		}
		else {
			rk[nsp - 1].len = lenk[nsp - 1];
			rk[nsp - 1].dat = malloc(sizeof(unsigned int)*rk[nsp - 1].len);
			for (int i = 0; i < rk[nsp - 1].len; i++) rk[nsp - 1].dat[i] = 1;
		}
	}
	/*! %synthesize r ------------------------------------------------------------------ */
	unsigned int * r_temp;
	for (int nsp = nsep; nsp >= 1; nsp--) {
		r_temp = separate_inv(rsep[nsp - 1].dat, (int *)rrw->dat, (int *)(rk[nsp - 1].dat), rsep[nsp - 1].len, thd[nsp]);
		free(rrw->dat);
		rrw->len = rsep[nsp - 1].len;
		rrw->dat = r_temp;
	}
#if 0
	内存释放：
		1 - thd
		2 - rsep 以及以下结构体成员
		3 - lenk
#endif // 0

		free(thd);
	for (int i = 0; i < nsep; i++) free(rsep[i].dat);
	free(rsep);
	free(lenk);

	return;
}
