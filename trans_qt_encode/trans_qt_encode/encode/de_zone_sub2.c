#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

DEC de_zone_sub2(uchar *biny, int lastBit,int lenbinbit)
{
	extern uchar *bin; //修复bin传不进de_runs_0seps中
	bin = biny;

	DEC ret;
	extern int ptr;
	int binBytePos = ptr & 7;


	DEC dec1 = decode_stationary_source_lenr(bin, 3, 1);
	DEC dec2 = decode_stationary_source_lenr(bin, 4, 1);
	int thd1 = (*(dec1.r));
	int thd2 = (*(dec2.r))*4;
	Uint8_Dat rSep2;
	de_runs_0seps(17, &rSep2);
	
	binBytePos = ptr & 7;
	int lB1 = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
	ptr++;
	int lenr1 = rSep2.len + 1 - lB1;
	//return ret;
	Uint32_Dat rs1;
	de_runs_1sep(lenr1, 4, &rs1, lenbinbit);

	GOLINV rsep1 = GolombInv(rs1.dat, lB1, rs1.len);
	DER rrw = de_r0(bin, sumBin(&rSep2), lenbinbit);
	Uint8_Dat rsep1Uint;
	rsep1Uint.dat = rsep1.z;
	rsep1Uint.len = rsep1.lenzbit;
	int lenk1 = rsep1.lenzbit - sumBin(&rsep1);
	int lenk2 = rSep2.len - sumBin(&rSep2);

	binBytePos = ptr & 7;
	int flg = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
	ptr++;

	DEC rk2,rk1;
	rk1.r = NULL;
	if (flg == 0)
	{
		rk2 = de_Kside(bin, lenk2, thd2, lenbinbit);
	}
	else
	{
		DES des = deSFcode(bin, 24);
		rk2 = decode_stationary_source_Nsym_lenr(bin, des.sym - 1, lenk2, thd2, lenbinbit);
	}

	if (thd1 > 1)
	{
		rk1 = de_Kside(bin, lenk1, thd1, lenbinbit);
	}
	else
	{
		rk1.lenr = lenk1;
		rk1.r = calloc(rk1.lenr, sizeof(unsigned int));
		for(int i = 0; i < rk1.lenr;i++)
		{ 
			rk1.r[i] = 1;
		}
	}
	Uint32_Dat rrw1;
	rrw1.len = rSep2.len;
	rrw1.dat =  separate_inv(rSep2.dat, rrw.r, rk2.r, rSep2.len, thd2);

	Uint32_Dat rrw2;
	rrw2.len = rsep1.lenzbit;
	rrw2.dat = separate_inv(rsep1.z, rrw1.dat,rk1.r, rrw2.len, thd1);
	GOLINV golRet = GolombInv(rrw2.dat, lastBit, rrw2.len);
	ret.z = golRet.z;
	ret.lenzbit = golRet.lenzbit;
	if (rSep2.dat != NULL)
	{
		free(rSep2.dat);
	}
	if (rs1.dat != NULL)
	{
		free(rs1.dat);
	}
	if (rsep1.z != NULL)
	{
		free(rsep1.z);
	}
	if (rrw.r != NULL)
	{
		free(rrw.r);
	}
	if (rk2.r != NULL)
	{
		free(rk2.r);
	}
	if (rk1.r != NULL)
	{
		free(rk1.r);
	}
	if (rrw1.dat != NULL)
	{
		free(rrw1.dat);
	}
	if (rrw2.dat != NULL)
	{
		free(rrw2.dat);
	}
	return ret;
}