#define  _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

void en_zone_sub2(Uint32_Dat*r)//biny得在调用此函数之前编入bin中去
{
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
	extern int ptr;
	extern uchar *bin;
	int binBytePos = ptr & 7;
	int* n = hist(r->dat, r->len); 
	int max = 0;
	for (int i = 0; i < r->len; i++)
	{
		if (max < r->dat[i])
		{
			max = r->dat[i];
		}
	}
	Int32_Dat nw,nk0, nk1;
	nw.dat = n;
	nw.len = max + 1;
	
	//float p[2] = {0.7031,0.8125};
	//ia = 1时
	int thd0 = find_thd(r->len, nw.dat, 0.7031);
	if (thd0 == 2)
	{
		if ((nw.dat[2] / (nw.dat[1]*1.0)) < 0.0945)
		{
			thd0 = 1;
		}
	}
	unsigned int*thd0Arr = calloc(1, sizeof(unsigned int));
	thd0Arr[0] = thd0;
	encode_stationary_source_bin(thd0Arr,1,3,0,0,0,0,0);
	//encode_stationary_source_bin(gol.r, gol.lenr, 1, 0, 0, 0, 0, 0);
	//ia = 1时
	SEP sep0 = separate0(r->dat, thd0, r->len, 1);
	nk0.dat = calloc(thd0, sizeof(unsigned int));
	nk0.len = thd0;

	Int32_Dat rk0;
	rk0.len = sep0.lrk;
	rk0.dat = calloc(rk0.len, sizeof(int));
	for (int i = 0; i < rk0.len; i++)
	{
		rk0.dat[i] = r->dat[i];
	}
	for (int i = 0; i < thd0; i++)
	{
		nk0.dat[i] = nw.dat[i + 1];
	}

	for (int i = thd0; i < nw.len; i++)
	{
		nw.dat[i - thd0] = nw.dat[i];
	}
	nw.len = nw.len - thd0;

	//ia = 2时
	int lenr = 0;
	for (int i = 1; i < nw.len; i++)
	{
		lenr += nw.dat[i];
	}
	int thd1 = find_thd(lenr, nw.dat, 0.8125);
	thd1 = round(thd1 /4.0) * 4;
	unsigned int*thd1Arr = calloc(1, sizeof(unsigned int));
	thd1Arr[0] = thd1 /4;
	encode_stationary_source_bin(thd1Arr, 1, 4, 0, 0, 0, 0, 0);

	SEP sep1 = separate0(sep0.rw, thd1, sep0.lrw, 1);
	nk1.dat = calloc(thd1, sizeof(unsigned int));
	nk1.len = thd1;
	for (int i = 0; i < thd1; i++)
	{
		nk1.dat[i] = nw.dat[i + 1];
	}
	for (int i = thd1; i < nw.len; i++)
	{
		nw.dat[i - thd1] = nw.dat[i];
	}
	nw.len = nw.len - thd1;

	Uint8_Dat sep21;
	sep21.dat = sep1.sep;
	sep21.len = sep1.lensepbit;

	en_runs_0seps(&sep21, 17);
	GOL gol = Golomb0(sep0.sep, sep0.lensepbit);

	binBytePos = ptr & 7;
	rem.a = gol.lastBit;//编lB
	rem.a = rem.a << (15 - binBytePos);
	bin[ptr >> 3] |= rem.b[1];
	bin[(ptr >> 3) + 1] |= rem.b[0];
	ptr++; binBytePos++;

	Uint32_Dat golR;
	golR.dat = gol.r;
	golR.len = gol.lenr;
	FILE *fp = fopen("golR.txt", "wb");
	fwrite(golR.dat, sizeof(unsigned int), golR.len, fp);
	fclose(fp);
	en_runs_1sep(&golR, 0.8125, 4, sep0.lensepbit, NULL);


	en_r0(sep1.rw, sep1.lrw);

	float prb1;
	prob_est(&nk1,&prb1);
	Uint32_Dat rk;
	rk.dat = sep0.rw;
	rk.len = sep1.lrk;
	int sumRk = 0;
	for (int i = 0; i < rk.len; i++)
	{
		sumRk += rk.dat[i];
	}
	float prb2 = rk.len / (1.0*(sumRk));
	float prob = 0;
	if (thd1 <= 6 && rk.len > 256)
	{
		prob = prb1;
	}
	else if (thd1 <= 12)
	{
		prob = (prb1 + prb2) / 2.0;
	}
	else
	{
		prob = prb2;
	}

	if (prob > 0.31)
	{
		ptr++; binBytePos++; //编0
		binBytePos &= 7;
		en_Kside(rk.dat, rk.len, thd1);
	}
	else
	{
		binBytePos = ptr & 7;
		rem.a = 1;//编1
		rem.a = rem.a << (15 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr++; binBytePos++;
		encode_stationary_source_Nsym(rk.dat, rk.len, prob, thd1, 1, 24);
	}
	if (thd0 > 1)
	{
		en_Kside(rk0.dat, rk0.len, thd0);
	}
	free(thd0Arr);
	free(thd1Arr);
	free(rk0.dat);
	if (n != NULL)
	{
		free(n);
	}
	if (nk0.dat != NULL)
	{
		free(nk0.dat);
	}
	if (nk1.dat != NULL)
	{
		free(nk1.dat);
	}
	if (sep0.sep != NULL)
	{
		free(sep0.sep);
	}
	if (sep0.rw != NULL)
	{
		free(sep0.rw);
	}
	if (sep1.sep != NULL)
	{
		free(sep1.sep);
	}
	if (sep1.rw != NULL)
	{
		free(sep1.rw);
	}
	if (gol.r != NULL)
	{
		free(gol.r);
	}
}