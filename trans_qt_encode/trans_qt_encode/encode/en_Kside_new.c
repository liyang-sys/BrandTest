#include <stdio.h>
#include "all.h"
#include "parameter_setting.h"
#include"encoding.h"

#define uint unsigned int
#define uchar unsigned char

//sumr可以用nk来算吗
//int
void en_Kside_new(int *rk, int l_rk, int *nk, int l_nk, int thd)
{
	extern int ptr;
	extern uchar *bin;
	int i, t0, t1, t;
	GOL g;
	ESSC e_cbook;
	SEP Separate0;
	int sumr = 0;
	float pk;
	int lens;                //中间值
	uchar codebook = 0;
	uchar *rA, *rB;    //把rk和rw转换为char的中间值
	Uint32_Dat nA;

	switch (thd) {
	case 2:
		rB = (uchar *)calloc((l_rk >> 3) + 1, sizeof(uchar));

		for (i = 0; i < l_rk; i++) {
			rk[i]--;
			t = i & 7;
			rB[i >> 3] |= rk[i] << (7 - t);
			sumr += rk[i];
		}
		/* pk=sum(rk-1)/length(rk); */
		pk = (float)(sumr) / l_rk;

		if (pk > 0.2421) {
			ptr++;   //把0移到bin
			en_KsideSub(rB, l_rk, sumr);
		}
		else {
			g = Golomb0(rB, l_rk);
			/* biny=[biny 1]; */
			t1 = ptr & 7;
			bin[ptr >> 3] |= (1 << (7 - t1));  //把1移到bin
			ptr++;


			e_cbook = encode_stationary_source_cbook(pk);

			/* codebook=uint8(abs( dec2bin(codebook-1,3) ) -48); */
			/* biny=[biny codebook]; */
			e_cbook.codebook--;
			for (i = 2; i >= 0; i--) {
				t0 = (e_cbook.codebook >> i) & 1;  //codebook的值直接赋给bin
				t1 = ptr & 7;
				bin[ptr >> 3] |= (t0 << (7 - t1));
				ptr++;
			}

			/* biny=[biny lastBit]; */
			t1 = ptr & 7;
			bin[ptr >> 3] |= (g.lastBit << (7 - t1));  //把lastBit移到bin
			ptr++;

			/* biny=[biny  bin]; */
			encode_stationary_source_bin(g.r, g.lenr, pk, e_cbook.k, e_cbook.m, e_cbook.m1, e_cbook.m2, e_cbook.cls);   //bin

			free(g.r);
		}
		free(rB);
		break;
	case 3:
		//printf("en_Kside 3:\n");///
		Separate0 = separate0(rk, 1, l_rk, 0);
		en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);


		rA = (uchar *)calloc((Separate0.lrw >> 3) + 1, sizeof(uchar));
		sumr = 0;
		for (i = 0; i < Separate0.lrw; i++) {
			Separate0.rw[i]--;
			t = i & 7;
			rA[i >> 3] |= Separate0.rw[i] << (7 - t);
			sumr += Separate0.rw[i];
		}

		en_KsideSub(rA, Separate0.lrw, sumr);

		free(Separate0.rw);
		free(Separate0.sep);
		free(rA);
		break;
	case 4:
		//printf("en_Kside 4:\n");///
		Separate0 = separate0(rk, 1, l_rk, 0);
		//enKsideSub =
		en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);
		lens = Separate0.lrw;

		Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
		//enKsideSub =
		en_KsideSub(Separate0.sep, lens, Separate0.lrw);

		rA = (uchar *)calloc((Separate0.lrw >> 3) + 1, sizeof(uchar));
		sumr = 0;
		for (i = 0; i < Separate0.lrw; i++) {
			Separate0.rw[i]--;
			t = i & 7;
			rA[i >> 3] |= Separate0.rw[i] << (7 - t);
			sumr += Separate0.rw[i];
		}

		//enKsideSub =
		en_KsideSub(rA, Separate0.lrw, sumr);

		free(Separate0.rw);
		free(Separate0.sep);
		free(rA);
		break;
	case 5:
		//printf("en_Kside 5:\n");///
		Separate0 = separate0(rk, 1, l_rk, 0);

		//enKsideSub =
		en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);
		lens = Separate0.lrw;

		Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
		//enKsideSub =
		en_KsideSub(Separate0.sep, lens, Separate0.lrw);

		lens = Separate0.lrw;
		Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
		//enKsideSub =
		en_KsideSub(Separate0.sep, lens, Separate0.lrw);

		rA = (uchar *)calloc((Separate0.lrw >> 3) + 1, sizeof(uchar));
		sumr = 0;
		for (i = 0; i < Separate0.lrw; i++)
		{
			Separate0.rw[i]--;
			t = i & 7;
			rA[i >> 3] |= Separate0.rw[i] << (7 - t);
			sumr += Separate0.rw[i];
		}

		//enKsideSub =
		en_KsideSub(rA, Separate0.lrw, sumr);

		free(Separate0.rw);
		free(Separate0.sep);
		free(rA);
		break;
	case 6:
		//printf("en_Kside 6:\n");///
		Separate0 = separate0(rk, 2, l_rk, 1);

		/*
		int s = 0;
		for (i = 0; i <= l_rk >> 3; i++)
		   s += Separate0.sep[i];*/

		   //enKsideSub =
		en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);
		lens = Separate0.lrw;

		rB = (uchar *)calloc((Separate0.lrk >> 3) + 1, sizeof(uchar));
		sumr = 0;
		for (i = 0; i < Separate0.lrk; i++) {
			rk[i]--;
			t = i & 7;
			rB[i >> 3] |= rk[i] << (7 - t);
			sumr += rk[i];
		}

		//enKsideSub =
		en_KsideSub(rB, Separate0.lrk, sumr);

		Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
		//enKsideSub =
		en_KsideSub(Separate0.sep, lens, Separate0.lrw);
		lens = Separate0.lrw;

		Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
		//enKsideSub =
		en_KsideSub(Separate0.sep, lens, Separate0.lrw);

		rA = (uchar *)calloc((Separate0.lrw >> 3) + 1, sizeof(uchar));
		sumr = 0;
		for (i = 0; i < Separate0.lrw; i++) {
			Separate0.rw[i]--;
			t = i & 7;
			rA[i >> 3] |= Separate0.rw[i] << (7 - t);
			sumr += Separate0.rw[i];
		}

		//enKsideSub =
		en_KsideSub(rA, Separate0.lrw, sumr);

		free(Separate0.rw);
		free(Separate0.sep);
		free(rA);
		break;
	default:
		//printf("en_Kside 7:\n");
		/* [sep,rA,rB]=separate(rk,2); rA=rA-2;  */
		Separate0 = separate0(rk, 2, l_rk, 1);
		/* nA=nk(3:thd); */
		nA.dat = (unsigned int *)(nk + 2);
		nA.len = thd - 3 + 1;
		/* [bin,len]=en_KsideSub(sep); biny=[biny bin]; */
		en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);  //sumr要重新求
		/* [bin,len]=en_KsideSub(rB-1); biny=[biny bin]; */
		rB = (uchar *)calloc((Separate0.lrk >> 3) + 1, sizeof(uchar));
		sumr = 0;
		for (i = 0; i < Separate0.lrk; i++) {
			rk[i]--;
			t = i & 7;
			rB[i >> 3] |= rk[i] << (7 - t);
			sumr += rk[i];
		}
		en_KsideSub(rB, Separate0.lrk, sumr);
		/* p=probility_estimation(rA,nA); */
		float p = 0;
		Uint32_Dat rA_temp;
		rA_temp.dat = Separate0.rw;
		rA_temp.len = Separate0.lrw;
		probility_estimation(&rA_temp, &nA, &p);
		//printf("ptr = %d\n", ptr);
		//printf("p =%f len rA = %d len nA = %d\n", p, rA_temp.len, nA.len);
		/* [bin,codebook]=encode_stationary_source_Nsym(rA,p,thd-2); */
		/* binck=encode_stationary_source(codebook+1, 4);  */
		/* biny = [biny binck bin];  */
		//printf("thd = %d\n", thd);
		//encode_stationary_source_Nsym(rA_temp.dat, rA_temp.len, p, thd - 2, 0, 4);
		encode_stationary_source_Nsym_in_en_Kside_new(rA_temp.dat, rA_temp.len, p, thd - 2, 4);

		free(Separate0.rw);
		free(Separate0.sep);
		free(rB);
		break;
	}
}
