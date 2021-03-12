#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

DE_S_SUB de_sub3d_sub2(uchar *bin, int len, int qctr, int lenbinbit)
{
	int thd_arr[24] = { 0,1, 2,3,4,6,8,12,	16,	24,	32,	48,	64,	96,	128,192,256,384,512	,768,1024,1536,2048	,3072 };
	DE_S_SUB ret;
	int Nidx = 23;
	extern int ptr;
	int binBytePos = ptr & 7;
	int Nsep = 0;
	Nsep = ((bin[ptr >> 3] >> (7 - binBytePos)) & 1)*4;
	ptr++;

	binBytePos = ptr & 7;
	Nsep = ((bin[ptr >> 3] >> (7 - binBytePos)) & 1)*2 + Nsep;
	ptr++;

	binBytePos = ptr & 7;
	Nsep = ((bin[ptr >> 3] >> (7 - binBytePos)) & 1) + Nsep + 1;
	ptr++;

	int* thd =(int*) calloc(Nsep + 1, sizeof(int));
	DE_S_S_SUB*sub_array = (DE_S_S_SUB*)calloc(Nsep + 1, sizeof(DE_S_S_SUB));

	for (int i = 1; i <= Nsep; i++)
	{
		DES thd_idx = deSFcode(bin, Nidx);
		thd[i] = thd_arr[thd_idx.sym];
		sub_array[i] = de_sub3d_sub2_sub(bin, thd[i], len, lenbinbit);

		Uint8_Dat temp;
		temp.dat = sub_array[i].z.dat;
		temp.len = sub_array[i].z.len;
		len = sumBin(&temp);
	}
	DER derN = de_r0(bin, len, lenbinbit);

	Int32_Dat cf;
	cf.dat = (int *)derN.r;
	cf.len = derN.lenr;

	ret.cfw.dat = (unsigned int*)calloc(derN.lenr, sizeof(unsigned int));
	ret.cfw.len = derN.lenr;
	for (int i = 0; i < derN.lenr; i++)
	{
		ret.cfw.dat[i] = derN.r[i];
	}

	for (int i = Nsep; i >= 1; i--)
	{
		Uint32_Dat temp;
		temp.dat = (unsigned int *)cf.dat;
		temp.len = cf.len;
		cf.dat = (int *)separate_inv(sub_array[i].z.dat, temp.dat, sub_array[i].cfk.dat, sub_array[i].z.len, thd[i]);
		cf.len = sub_array[i].z.len;
		free(temp.dat);
	}

	if (qctr == 1)
	{
		for (int i = 0; i < cf.len; i++)
		{
			(cf.dat[i])--;
		}
		int lensn = 0;
		for (int i = 0; i < cf.len; i++)
		{
			if (cf.dat[i] < 0 || cf.dat[i] > 0)
			{
				lensn++;
			}
		}
		Int32_Dat sn;
		sn.dat = calloc(lensn, sizeof(int));
		sn.len = lensn;

		for (int i = ptr; i <= (ptr + lensn -1); i++)
		{
			binBytePos = i & 7;
			if ((bin[i >> 3] >> (7 - binBytePos)) & 1)
			{
				sn.dat[i - ptr] = 1;
			}
			else
			{
				sn.dat[i - ptr] = -1;
			}
		}
		ptr += lensn;
		int ns = 0;
		for (int i = 0; i < cf.len; i++)
		{
			if (cf.dat[i] != 0)
			{
				cf.dat[i] = (cf.dat[i]) * (sn.dat[ns]);
				ns++;
			}
		}

		ret.cf.dat = cf.dat;
		ret.cf.len = cf.len;
		ret.sn.dat = sn.dat;
		ret.sn.len = sn.len;

	}
	else
	{
		int lensn = cf.len;
		Int32_Dat sn;
		sn.dat = calloc(lensn, sizeof(int));
		sn.len = lensn;

		for (int i = ptr; i <= (ptr + lensn - 1); i++)
		{
			binBytePos = i & 7;
			if ((bin[i >> 3] >> (7 - binBytePos)) & 1)
			{
				sn.dat[i - ptr] = 1;
			}
			else
			{
				sn.dat[i - ptr] = -1;
			}
		}
		ptr += lensn;
		for (int i = 0; i < cf.len; i++)
		{

			cf.dat[i] = (cf.dat[i]) * (sn.dat[i]);
			
		}
		ret.cf.dat = cf.dat;
		ret.cf.len = cf.len;
		ret.sn.dat = sn.dat;
		ret.sn.len = sn.len;
	}

	free(thd);
	for (int i = 1; i <= Nsep; i++)
	{
		free(sub_array[i].cfk.dat);
		free(sub_array[i].z.dat);
	}
	return ret;
}