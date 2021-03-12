#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

DEC de_zone_new(uchar *biny,int lenbinbit)
{
	uchar* bin;
	bin = biny;
	DEC ret;
	extern int ptr;
	int binBytePos = ptr & 7;
	int th[25] = { 0,1,2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024,1536,2048,3072,4096 };

	int lB = 0;
	lB = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
	ptr++;

	DEC idx = decode_stationary_source_lenr(bin, 3, 1);

	int thd1 = th[idx.r[0]];
	int lBs = 0;
	binBytePos = ptr & 7;
	lBs = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
	ptr++;

	int flg = 0;
	binBytePos = ptr & 7;
	flg = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
	ptr++;
	Uint32_Dat rs;
	if (flg == 1)
	{
		Uint32_Dat rsNoneed;
		rsNoneed.dat = NULL;
		de_zone_sub_new(bin, lenbinbit, &rs, &rsNoneed);
		if (rsNoneed.dat != NULL)
		{
			free(rsNoneed.dat);
		}
	}
	else
	{
		de_runs_2sepsA(17, &rs, lenbinbit);
	}

	GOLINV sep1 = GolombInv(rs.dat, lBs, rs.len);
	int minsep1 = 0 , maxsep1 = 0;
	for (int i = 0; i < ((sep1.lenzbit / 8) + 1); i++)
	{
		if (sep1.z[i] > 0)
		{
			maxsep1 = 1;
			break;
		}
	}
	Uint8_Dat binZ;
	binZ.dat = sep1.z;
	binZ.len = sep1.lenzbit;
	int lenk1 = sep1.lenzbit - sumBin(&binZ);
	Uint32_Dat rk1;

	if (minsep1 == 0)
	{
		if (thd1 > 1)
		{
			binBytePos = ptr & 7;
			flg = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
			ptr++;
			if (flg > 0)
			{
				int codebook = 0;
				for (int i = 0; i < 4; i++)
				{
					codebook = codebook << 1;
					binBytePos = ptr & 7;
					//codebook = (bin[ptr >> 3] >> (7 - binBytePos)) & 1 + codebook;
					codebook = ((bin[ptr >> 3] >> (7 - binBytePos)) & 1) + codebook;
					ptr++;

				}
				DEC rkDec = decode_stationary_source_Nsym_lenr(bin, codebook, lenk1, thd1, lenbinbit);
				rk1.dat = rkDec.r;
				rk1.len = rkDec.lenr;
			}
			else
			{
				DEC rkDec = de_Kside_new(bin, lenk1, thd1, lenbinbit);
				rk1.dat = rkDec.r;
				rk1.len = rkDec.lenr;
			}
		}
		else
		{
			rk1.dat = calloc(lenk1, sizeof(unsigned int));
			rk1.len = lenk1;
			for (int i = 0; i < rk1.len; i++)
			{
				rk1.dat[i] = 1;
			}
		}
	}
	Uint32_Dat rs2,rkNoNeed,rk2,rw2, rk2NoNeed,rw1, rk3NoNeed;
	rkNoNeed.dat = NULL;
	rk2NoNeed.dat = NULL;
	rk3NoNeed.dat = NULL;
	rs2.dat = NULL;
	rk2.dat = NULL;
	rw1.dat = NULL;
	rw2.dat = NULL;

	if (maxsep1 == 1)
	{
		binBytePos = ptr & 7;
		flg = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;

		if (flg == 1)
		{
			DEC idx2 = decode_stationary_source_lenr(bin, 3, 1);
			int thd2 = th[idx2.r[0]];
			binBytePos = ptr & 7;
			lBs = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
			ptr++;
			de_zone_sub_new(bin, lenbinbit, &rs2, &rkNoNeed);

			GOLINV sep2 = GolombInv(rs2.dat, lBs, rs2.len);
			Uint8_Dat binZ2;
			binZ2.dat = sep2.z;
			binZ2.len = sep2.lenzbit;

			int lenk2 = sep2.lenzbit - sumBin(&binZ2);

			if (thd2 > 1)
			{
				binBytePos = ptr & 7;
				flg = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
				ptr++;
				if (flg > 0)
				{
					int codebook = 0;
					for (int i = 0; i < 4; i++)
					{
						codebook = codebook << 1;
						binBytePos = ptr & 7;
						//codebook = (bin[ptr >> 3] >> (7 - binBytePos)) & 1 + codebook;
						codebook = ((bin[ptr >> 3] >> (7 - binBytePos)) & 1) + codebook;
						ptr++;

					}
					DEC rkDec = decode_stationary_source_Nsym_lenr(bin, codebook, lenk2, thd2, lenbinbit);
					rk2.dat = rkDec.r;
					rk2.len = rkDec.lenr;
				}
				else
				{
					DEC rkDec = de_Kside_new(bin, lenk2, thd2, lenbinbit);
					rk2.dat = rkDec.r;
					rk2.len = rkDec.lenr;
				}
			}
			else
			{
				rk2.dat = calloc(lenk2, sizeof(unsigned int));
				rk2.len = lenk2;
				for (int i = 0; i < rk2.len; i++)
				{
					rk2.dat[i] = 1;
				}
			}
			de_zone_sub_new(bin, lenbinbit, &rw2, &rk2NoNeed);
			rw1.dat = separate_inv(sep2.z, rw2.dat, rk2.dat, sep2.lenzbit, thd2);
			rw1.len = sep2.lenzbit;

			if (idx2.r != NULL)
			{
				free(idx2.r);
			}
			if (sep2.z != NULL)
			{
				free(sep2.z);
			}
		}
		else
		{
			de_zone_sub_new(bin, lenbinbit, &rw1, &rk3NoNeed);
		}
	}

	Uint32_Dat r;
	if (minsep1 == 0)
	{
		if (maxsep1 == 1)
		{
			r.dat = separate_inv(sep1.z, rw1.dat, rk1.dat,sep1.lenzbit, thd1);
			r.len = sep1.lenzbit;
		}
		else
		{
			r.dat = rk1.dat;
			r.len = rk1.len;
		}
	}
	else
	{
		for (int i = 0; i < rw1.len; i++)
		{
			(rw1.dat[i]) += thd1;
		}
		r.dat = rw1.dat;
		r.len = rw1.len;
	}
	GOLINV gol = GolombInv(r.dat, lB, r.len);
	ret.z = gol.z;
	ret.lenzbit = gol.lenzbit;

	if (idx.r != NULL)
	{
		free(idx.r);
	}
	if (rs.dat != NULL)
	{
		free(rs.dat);
	}
	if (rk1.dat != NULL)
	{
		free(rk1.dat);
	}
	if (rkNoNeed.dat != NULL)
	{
		free(rkNoNeed.dat);
	}
	if (rs2.dat != NULL)
	{
		free(rs2.dat);
	}
	if (rk2.dat != NULL)
	{
		free(rk2.dat);
	}
	if (rk2NoNeed.dat != NULL)
	{
		free(rk2NoNeed.dat);
	}
	if (rw2.dat != NULL)
	{
		free(rw2.dat);
	}
	if (rw1.dat != NULL)
	{
		free(rw1.dat);
	}
	if (rk3NoNeed.dat != NULL)
	{
		free(rk2NoNeed.dat);
	}
	if (r.dat != NULL)
	{
		if ((r.dat != rk1.dat) && (r.dat != rw1.dat))
		{
			free(r.dat);
		}
	}
	return ret;
}