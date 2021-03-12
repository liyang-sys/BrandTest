#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

int snbin_add(uchar *bin1, uchar *snbin, int lenbin2bit)
{
	union data {
		uint a;
		uchar b[4];
	} rem;
	extern int ptr;
	int i, x, p;
	int lenbin2;
	x = ptr & 7;
	p = ptr;

	i = 0;
	lenbin2 = len_bit_to_byte(lenbin2bit);
	while (i < lenbin2)
	{
		rem.b[3] = (uchar)(~snbin[i]);
		rem.b[2] = (uchar)(~snbin[i + 1]);
		rem.b[1] = (uchar)(~snbin[i + 2]);
		rem.b[0] = (uchar)(~snbin[i + 3]);

		rem.a = rem.a >> x;
		bin1[p >> 3] |= rem.b[3];
		bin1[(p >> 3) + 1] = rem.b[2];
		bin1[(p >> 3) + 2] = rem.b[1];
		bin1[(p >> 3) + 3] = rem.b[0];

		i += 3;
		p += 24;
	}
	ptr += lenbin2bit;
	//x = ptr & 7;
	//bin1[ptr >> 3] = (bin1[ptr >> 3] >> (7 - x)) << (7 - x); 
	//bin1[(ptr >> 3) + 1] = 0;
	//bin1[(ptr >> 3) + 2] = 0;
	//bin1[(ptr >> 3) + 3] = 0;

	//对多出来的bin置0
	int lastidx = ptr >> 3;	
	int resbin = ptr - (lastidx << 3);
	bin1[lastidx] >>= (8 - resbin);
	bin1[lastidx] <<= (8 - resbin);

	return 0;
}


void en_sub3d_sub2(float* cf0, Uint8_Dat* sn, int lencf0, int maxcf0)
{

	/* Hint:nc以及thd的0索引保留不使用 */

	extern unsigned char *bin;
	extern int ptr;
	int tempPtr = ptr;
	ptr += 3;

	union data {
		uint a;
		uchar b[4];
	} rem;
	int x = 0;
	SFC se;

	int maxcf = 0;
	//unsigned int len = 0;
	unsigned int len = lencf0;//修改的
	unsigned int *nc0 = NULL;
	unsigned int Nidx = 0;
	float p = 0, p1 = 0, p2 = 0.8;
	unsigned int lenw = 0, lenc = 0;
	Int32_Dat cf;
	int iter = 0, idx0 = 0, crc = 20;
	int thd1 = 0;

	int thdidx = 0, lB = 0;
	Uint32_Dat r, nr, cf1, cf1k, nc1, nc1k;
	float pr = 0, scr = 0, scrc = 0;
	Uint8_Dat z;
	Int32_Dat thd_arr_temp;
	Uint32_Dat nc_temp;
	int nc_len = 0;

	cf.dat = (int*)malloc(lencf0*sizeof(int));

	/*!  %<=== treat cf0 as runs! */
	/* maxcf=double(maxcf0+1); */
	maxcf = maxcf0 + 1;
	nc0 = (unsigned int *)malloc(sizeof(unsigned int)*(maxcf + 1));
	if (nc0 == NULL) {
		perror("n");
		exit(-1);
	}
	memset(nc0, 0, sizeof(unsigned int)*(maxcf + 1));

	/* cf0=cf0+1; nc=hist(double(cf0),[1:maxcf]); */
	for (int i = 0; i < lencf0; i++) {
		cf.dat[i] = (int)cf0[i] + 1;
		nc0[(unsigned int)cf.dat[i]]++;
	}

	int thd_arr[24] = { 0, 1, 2, 3, 4, 6, 8,
					  12, 16, 24, 32, 48, 64, 96, 128,
					  192, 256, 384, 512, 768, 1024, 1536, 2048, 3072 };
	Nidx = 23;
	p = 0.7; p1 = 0.7; p2 = 0.8;
	p = p;   p1 = p1; p2 = p2;   // 保留，防止编译器报警告

	/*! %first round ========================================================================== */
	lenw = len;
	lenc = 1600;
	iter = 1;
	cf.len = lencf0;
	idx0 = 1;
	maxcf = maxcf;
	crc = 20;

	nc_len = maxcf;
	unsigned int *temp = (unsigned int *)malloc(sizeof(unsigned int)*maxcf);
	memcpy(temp, nc0 + 1, sizeof(unsigned int)*maxcf);
	free(nc0);
	nc0 = temp;

	while ((lenw > lenc || crc >= 12) && iter <= 8 && maxcf >= 2) {
		thd_arr_temp.dat = thd_arr + 1;
		thd_arr_temp.len = Nidx;
		nc_temp.dat = nc0;
		nc_temp.len = nc_len;
		find_thd4subs(&cf, &thd_arr_temp, &nc_temp, idx0, &thdidx, &r, &lB, &pr, &nr, &scr, &scrc, &cf1, &cf1k, &nc1, &nc1k, &z);
		thd1 = thd_arr[thdidx];
		thd1 = thd1;                 // 保留，防止编译器报警告
		/* idx0=thdidx; */
		idx0 = thdidx;
		/* biny=[biny SFcode(thdidx,Nidx)]; */
		x = ptr & 7;
		se = SFcode(thdidx, Nidx);
		rem.a = se.code;
		rem.a = rem.a << (16 - x - se.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += se.lb; x += se.lb; x &= 7;

		Uint32_Dat nrPlus;
		nrPlus.dat = nr.dat + 1;
		nrPlus.len = nr.len;

		en_sub3d_sub2_sub(&z, &r, lB, pr, &nrPlus, scr, scrc, &cf1k, &nc1k, thd1);

		/* iter=iter+1; */
		iter += 1;
		/* lenw=length(cf1); */
		lenw = cf1.len;
		/* cf=cf1; */
		free(cf.dat);
		cf.dat = (int *)cf1.dat;
		cf.len = cf1.len;

		/* nc0=nc1; */
		free(nc0);
		nc0 = nc1.dat;
		nc_len = nc1.len;
		maxcf = 0;
		for (int i = 0; i < cf.len; i++) maxcf = maxcf > cf.dat[i] ? maxcf : cf.dat[i];

		//代填
		if (maxcf >= 2) {
			free(nc0);//修改
			k_criterion((void *)&cf, NULL, &crc, NULL, NULL, NULL, &nc_temp, NULL, 0);
			nc0 = nc_temp.dat;
			nc_len = nc_temp.len;
			for (int i = 0; i < nc_len; i++)
			{
				nc0[i] = nc0[i + 1];
			}
			//free(nc_temp.dat);//修改
			//printf("crc = %d\n", crc);
		}
		else {
			break;
		}
	}
	iter += -1;

	rem.a = iter - 1;
	rem.a = rem.a << (16 - (tempPtr & 7) - 3);
	bin[tempPtr >> 3] |= rem.b[1];
	bin[(tempPtr >> 3) + 1] |= rem.b[0];

	en_r0((unsigned int *)cf.dat, cf.len);
	x &= ptr;

	snbin_add(bin,sn->dat , sn->len);
	free(nc0);

	return;
}
