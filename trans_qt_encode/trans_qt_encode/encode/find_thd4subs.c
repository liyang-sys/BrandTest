#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

/**
  * @brief
  * @param  cf0[IN]
  * @param  thd_arr[IN] : 从索引0开始使用
  * @param  nc[IN]      : 从索引0开始使用
  * @param  idx0[IN]
  * @param  idx[OUT]
  * @param  r[OUT]
  * @param  lB[OUT]
  * @param  pr[OUT]
  * @param  nr[OUT] 数据得从索引1开始使用
  * @param  sumcr[OUT]
  * @param  sumcrc[OUT]
  * @param  cf1k[OUT]
  * @param  nc1[OUT]
  * @param  nc1k[OUT]
  * @param  z1[OUT]      : 存的是二进制流的长度
  * @retval none
  * @note
  */
void find_thd4subs(Int32_Dat* cf0, Int32_Dat* thd_arr, Uint32_Dat* nc, int idx0, int *idx, Uint32_Dat* r, int *lB,
	float *pr, Uint32_Dat* nr, float* sumcr, float* sumcrc, Uint32_Dat* cf1, Uint32_Dat* cf1k,
	Uint32_Dat* nc1, Uint32_Dat* nc1k, Uint8_Dat* z1)
{
	float p = 0.5;
	*idx = 0;
	int   thd = 0;
	float cr = 0;
	int   crc = 0;
	float cri = 0;
	int   crci = 0;
	int   lBi = 0;
	Uint32_Dat ri;
	float maxcr = 0;
	float sumcr0 = 0;
	float sumcrc0 = 0;
	int   flg = 0;
	int   max_idx = 0;
	int   sum_nc1k = 0;
	int   sum_nc = 0;
	SEP   sep_temp;
	Int32_Dat cf0_bak;
	maxcr = maxcr;

	p = 0.5;
	find_the_idx_round(nc, p, thd_arr->dat, idx);
	if (idx0 > *idx && thd_arr->dat[*idx] < nc->len) {
		(*idx)++;
	}

	thd = thd_arr->dat[(unsigned int)((*idx) - 1)];
	cf0_bak.len = cf0->len;
	cf0_bak.dat = (int *)malloc(sizeof(int)*(cf0->len));
	memcpy(cf0_bak.dat, cf0->dat, sizeof(int)*(cf0->len));
	sep_temp = separate0(cf0_bak.dat, thd, cf0_bak.len, 1);
	z1->dat = sep_temp.sep;
	z1->len = sep_temp.lensepbit;
	cf1->dat = sep_temp.rw;
	cf1->len = sep_temp.lrw;
	cf1k->dat = (unsigned int*)cf0_bak.dat;
	cf0_bak.dat = NULL;
	cf1k->len = sep_temp.lrk;

	/* nc1k=nc(1:thd); */
	nc1k->len = thd;
	nc1k->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1k->len);
	memcpy(nc1k->dat, nc->dat, sizeof(unsigned int)*nc1k->len);

	/* nc1=nc(thd+1:length(nc)); */
	nc1->len = nc->len - thd;
	nc1->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1->len);
	memcpy(nc1->dat, &(nc->dat[thd]), sizeof(unsigned int)*nc1->len);

	/* [cr,crc,r,lB,pr,nr]=k_criterion(z1); */
	k_criterion(z1, &cr, &crc, r, lB, pr, nr, NULL, 1);

	/* [cri,crci,ri,lBi]=k_criterion(1-z1); */
	for (int i = 0; i < (z1->len + 7) / 8; i++) {
		z1->dat[i] = ~z1->dat[i];
	}
	k_criterion(z1, &cri, &crci, &ri, &lBi, NULL, NULL, NULL, 1);

	/* maxcr=max([cr cri]); */
	maxcr = cr > cri ? cr : cri;

	sumcr0 = 0;
	*sumcr = cr + cri;
	*sumcrc = crc + crci;

	flg = 0;
	/* max_idx=numel(thd_arr); */
	max_idx = thd_arr->len;

	/* sum(nc) */
	sum_nc = 0;
	for (int i = 0; i < nc->len; i++)
		sum_nc += nc->dat[i];
	/* sum(nc1k) */
	sum_nc1k = 0;
	for (int i = 0; i < nc1k->len; i++)
		sum_nc1k += nc1k->dat[i];

	while (*sumcr > sumcr0 && (float)sum_nc1k / sum_nc < 0.86 && *idx < max_idx && thd_arr->dat[*idx] < nc->len) {
		flg = 1;
		sumcr0 = *sumcr;
		/* 好像这句话没有什么作用，但还是保留下来吧 */
		sumcrc0 = *sumcrc;
		sumcrc0 = sumcrc0;
		*idx += 1;
		thd = thd_arr->dat[(unsigned int)(*idx - 1)];

		/* [z1,cf1,cf1k]=separate(cf0,thd); */
		/* cf1=cf1-thd;  */
		free(z1->dat);
		free(cf1->dat);
		free(cf1k->dat);
		cf0_bak.len = cf0->len;
		cf0_bak.dat = (int *)malloc(sizeof(int)*(cf0->len));
		memcpy(cf0_bak.dat, cf0->dat, sizeof(int)*(cf0->len));
		sep_temp = separate0(cf0_bak.dat, thd, cf0_bak.len, 1);
		z1->dat = sep_temp.sep;
		z1->len = sep_temp.lensepbit;
		cf1->dat = sep_temp.rw;
		cf1->len = sep_temp.lrw;
		cf1k->dat = (unsigned int*)cf0_bak.dat;
		cf0_bak.dat = NULL;
		cf1k->len = sep_temp.lrk;

		/* nc1k=nc(1:thd); */
		nc1k->len = thd;
		free(nc1k->dat);
		nc1k->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1k->len);
		memcpy(nc1k->dat, nc->dat, sizeof(unsigned int)*nc1k->len);

		sum_nc1k = 0;//缺少，由2021-03-06添加，zw。
		for (int i = 0; i < nc1k->len; i++)
			sum_nc1k += nc1k->dat[i];

		/* nc1=nc(thd+1:length(nc)); */
		nc1->len = nc->len - thd;
		free(nc1->dat);
		nc1->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1->len);
		memcpy(nc1->dat, &(nc->dat[thd]), sizeof(unsigned int)*nc1->len);

		/* [cr,crc,r,lB,pr,nr]=k_criterion(z1); */
		free(r->dat);
		free(nr->dat);
		k_criterion(z1, &cr, &crc, r, lB, pr, nr, NULL, 1);

		/* [cri,crci,ri,lBi]=k_criterion(1-z1); */
		free(ri.dat);
		for (int i = 0; i < (z1->len + 7) / 8; i++) {
			z1->dat[i] = ~z1->dat[i];
		}
		k_criterion(z1, &cri, &crci, &ri, &lBi, NULL, NULL, NULL, 1);

		/* maxcr=max([cr cri]); */
		maxcr = cr > cri ? cr : cri;

		*sumcr = cr + cri;
		*sumcrc = crc + crci;
	}

	if (flg == 1) {
		*idx += -1;
		thd = thd_arr->dat[(unsigned int)(*idx - 1)];
		/* [z1,cf1,cf1k]=separate(cf0,thd); */
		/* cf1=cf1-thd;  */
		free(z1->dat);
		free(cf1->dat);
		free(cf1k->dat);
		cf0_bak.len = cf0->len;
		cf0_bak.dat = (int *)malloc(sizeof(int)*(cf0->len));
		memcpy(cf0_bak.dat, cf0->dat, sizeof(int)*(cf0->len));
		sep_temp = separate0(cf0_bak.dat, thd, cf0_bak.len, 1);
		z1->dat = sep_temp.sep;
		z1->len = sep_temp.lensepbit;
		cf1->dat = sep_temp.rw;
		cf1->len = sep_temp.lrw;
		cf1k->dat = (unsigned int*)cf0_bak.dat;
		cf0_bak.dat = NULL;
		cf1k->len = sep_temp.lrk;

		/* nc1k=nc(1:thd); */
		nc1k->len = thd;
		free(nc1k->dat);
		nc1k->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1k->len);
		memcpy(nc1k->dat, nc->dat, sizeof(unsigned int)*nc1k->len);
		/* nc1=nc(thd+1:length(nc)); */
		nc1->len = nc->len - thd;
		free(nc1->dat);
		nc1->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1->len);
		memcpy(nc1->dat, &(nc->dat[thd]), sizeof(unsigned int)*nc1->len);

		/* [cr,crc,r,lB,pr,nr]=k_criterion(z1); */
		free(r->dat);
		free(nr->dat);
		k_criterion(z1, &cr, &crc, r, lB, pr, nr, NULL, 1);

		/* [cri,crci,ri,lBi]=k_criterion(1-z1); */
		free(ri.dat);
		for (int i = 0; i < (z1->len + 7) / 8; i++) {
			z1->dat[i] = ~z1->dat[i];
		}
		k_criterion(z1, &cri, &crci, &ri, &lBi, NULL, NULL, NULL, 1);

		/* maxcr=max([cr cri]); */
		maxcr = cr > cri ? cr : cri;

		*sumcr = cr + cri;
		*sumcrc = crc + crci;
	}

	free(ri.dat);
	ri.dat = NULL;
	//添加的
	for (int i = 0; i < (z1->len + 7) / 8; i++) {
		z1->dat[i] = ~z1->dat[i];
	}
	return;
}
