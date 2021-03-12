#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

void en_zone_new(Uint32_Dat*r, Uint32_Dat*nc)
{
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
	extern int ptr;
	extern uchar *bin;
	int binBytePos = ptr & 7;

	int th[25] = {0,1,2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024,1536,2048,3072,4096 };
	float p1 = 0.7, p2 = 0.8;
	int Ln = nc->len;
	
	int idx = find_thd_idx(nc,p1,th);//问一下find_thd_idx区别
	int thd1 = th[idx];
	
	Uint32_Dat idxU32;
	idxU32.len = 1;
	idxU32.dat = calloc(1, sizeof(unsigned int));
	idxU32.dat[0] = idx;
	encode_stationary_source_bin(idxU32.dat, idxU32.len, 3, 0, 0, 0, 0,0);

	SEP sep1 = separate0(r->dat, thd1, r->len, 1);
	Uint32_Dat nck;
	nck.dat = calloc(thd1, sizeof(unsigned int));
	nck.len = thd1;
	for (int i = 0; i < thd1; i++)
	{
		nck.dat[i] = nc->dat[i];
	}

	Uint32_Dat ncw;
	ncw.len = Ln-thd1;
	ncw.dat = calloc(ncw.len, sizeof(unsigned int));
	for (int i = thd1; i < Ln; i++)
	{
		ncw.dat[i-thd1] = nc->dat[i];
	}
	Ln -= thd1;

	float crs = 0;
	Uint32_Dat rs;
	rs.dat = NULL;
	int lBs = 0;
	Uint8_Dat z;
	z.dat = sep1.sep;
	z.len = sep1.lensepbit;
	k_criterion(&z, &crs, NULL, &rs, &lBs, NULL, NULL, NULL, 1);
	if (crs > 60)
	{
		binBytePos =ptr & 7;
		rem.a = lBs;//编lBs
		rem.a = rem.a << (15 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr++; binBytePos++;

		binBytePos = ptr & 7;
		rem.a = 1;//编1
		rem.a = rem.a << (15 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr++; binBytePos++;

		SEP* sepTemp;
		en_zone_sub_new(&rs, crs, &sepTemp);
		//en_zone_sub_new
		if (sepTemp->nk != NULL)
		{
			free(sepTemp->nk);
		}
		if (sepTemp->rk != NULL)
		{
			free(sepTemp->rk);
		}
		if (sepTemp->rw != NULL)
		{
			free(sepTemp->rw);
		}
		if (sepTemp->sep != NULL)
		{
			free(sepTemp->sep);
		}
	}
	else
	{
		binBytePos = ptr & 7;
		rem.a = lBs;//编lBs
		rem.a = rem.a << (15 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr++; binBytePos++;
		ptr++;//编0
		en_runs_2sepsA(&rs, p2, 17, NULL);
		//en_runs_2sepsA
	}

	if (thd1 > 1 && (r->len != 0))
	{
		Int32_Dat rk1;
		rk1.dat = r->dat;
		rk1.len = sep1.lrk;

		Int32_Dat nckInt32;
		nckInt32.dat = nck.dat;
		nckInt32.len = nck.len;

		float prob;
		probility_estimation(&rk1, &nckInt32, &prob);
		if ((prob > 0.15) && (thd1 <= 16))
		{
			ptr++;//编0
			en_Kside_new(rk1.dat, rk1.len, nckInt32.dat, nckInt32.len, thd1);

		}
		else
		{
			binBytePos = ptr & 7;
			rem.a = 1;//编1
			rem.a = rem.a << (15 - binBytePos);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr++; binBytePos++;

			encode_stationary_source_Nsym(rk1.dat, rk1.len, prob, thd1, 0,0);//直接 0 0 传参
		}
	}
	Uint32_Dat rw;
	rw.dat = sep1.rw;
	rw.len = sep1.lrw;
	Uint32_Dat rwCopy;
	rwCopy.dat = calloc(rw.len, sizeof(unsigned int));
	if (rwCopy.dat == NULL)
	{
		printf("rwCopy分配内存失败\n");
	}
	rwCopy.len = rw.len;
	for (int i = 0; i < rwCopy.len; i++)
	{
		rwCopy.dat[i] = rw.dat[i];
	}
	if (rw.len != 0)
	{
		idx = find_thd_idx(&ncw, p1, th);
		int thd2 = th[idx];
		SEP sep2 = separate0(rw.dat, thd2, rw.len, 1);
		Uint32_Dat nck2;
		nck2.dat = calloc(thd2, sizeof(unsigned int));
		nck2.len = thd2;
		for (int i = 0; i < thd2; i++)
		{
			nck2.dat[i] = ncw.dat[i];
		}

		float crs2 = 0;
		Uint32_Dat rs2;
		rs2.dat = NULL;
		int lBs2 = 0;
		Uint8_Dat z2;
		z2.dat = sep2.sep;
		z2.len = sep2.lensepbit;
		k_criterion(&z2, &crs2, NULL, &rs2, &lBs2, NULL, NULL, NULL, 1);
		if (crs2 > 50)
		{
			binBytePos = ptr & 7;
			rem.a = 1;//编1
			rem.a = rem.a << (15 - binBytePos);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr++; binBytePos++;

			Uint32_Dat idx2U32;
			idx2U32.len = 1;
			idx2U32.dat = calloc(1, sizeof(unsigned int));
			idx2U32.dat[0] = idx;
			
			encode_stationary_source_bin(idx2U32.dat, idx2U32.len, 3, 0, 0, 0, 0, 0);//编bin_idx2
			
			binBytePos = ptr & 7; 
			rem.a = lBs2;//编lBs 
			rem.a = rem.a << (15 - binBytePos);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr++; binBytePos++;

			SEP* sepTemp;
			en_zone_sub_new(&rs2, crs2, &sepTemp);//编bins2
			//[bins2]=en_zone_sub_new(rs,crs,lBs);
			if (sepTemp->nk != NULL)
			{
				free(sepTemp->nk);
			}
			if (sepTemp->rk != NULL)
			{
				free(sepTemp->rk);
			}
			if (sepTemp->rw != NULL)
			{
				free(sepTemp->rw);
			}
			if (sepTemp->sep != NULL)
			{
				free(sepTemp->sep);
			}

			if (thd2 > 1)
			{
				Uint32_Dat rk2;
				rk2.dat = rw.dat;
				rk2.len = sep2.lrk;
				float prob2;
				probility_estimation(&rk2, &nck2, &prob2);
				if ((prob2 > 0.15) & (thd2 <= 16))
				{
					ptr++;//编0
					en_Kside_new(rk2.dat, rk2.len, nck2.dat, nck2.len, thd2);
				}
				else
				{
					binBytePos  = ptr& 7;
					rem.a = 1;//编1
					rem.a = rem.a << (15 - binBytePos);
					bin[ptr >> 3] |= rem.b[1];
					bin[(ptr >> 3) + 1] |= rem.b[0];
					ptr++; binBytePos++;
					encode_stationary_source_Nsym(rk2.dat, rk2.len, prob2, thd2, 0, 0);//直接 0 0 传参
				}
			}
			Uint32_Dat rw2;
			rw2.dat = sep2.rw;
			rw2.len = sep2.lrw;
		

			float crw2;
			k_criterion(&rw2, &crw2, NULL, NULL, NULL, NULL, NULL, NULL, 0);
			SEP* sepTemp1;
			en_zone_sub_new(&rw2, crw2, &sepTemp1);

			//K_Criterion()
			if (sepTemp1->nk != NULL)
			{
				free(sepTemp1->nk);
			}
			if (sepTemp1->rk != NULL)
			{
				free(sepTemp1->rk);
			}
			if (sepTemp1->rw != NULL)
			{
				free(sepTemp1->rw);
			}
			if (sepTemp1->sep != NULL)
			{
				free(sepTemp1->sep);
			}

			if (idx2U32.dat != NULL)
			{
				free(idx2U32.dat);
			}
			//if (rw2Copy.dat != NULL)
			//{
			//	free(rw2Copy.dat);
			//}
		}
		else
		{
			ptr++;//编0
			float crw;
			k_criterion(&rwCopy, &crw, NULL,NULL, NULL, NULL, NULL, NULL, 0);
			SEP* sepTemp;

			FILE* fp = fopen("rwCopy.txt", "wb");
			fwrite(rwCopy.dat, sizeof(unsigned int), rwCopy.len, fp);
			fclose(fp);
			en_zone_sub_new(&rwCopy, crw, &sepTemp);
			
			//K_Criterion()
			if (sepTemp->nk != NULL)
			{
				free(sepTemp->nk);
			}
			if (sepTemp->rk != NULL)
			{
				free(sepTemp->rk);
			}
			if (sepTemp->rw != NULL)
			{
				free(sepTemp->rw);
			}
			if (sepTemp->sep != NULL)
			{
				free(sepTemp->sep);
			}
		}
		if (nck2.dat != NULL)
		{
			free(nck2.dat);
		}
		//if (sep2.sep != NULL)
		//{
		//	free(sep2.sep);
		//}
		//if (sep2.rw != NULL)
		//{
		//	free(sep2.rw);
		//}
		if (rs2.dat != NULL)
		{
			free(rs2.dat);
		}
		if (rwCopy.dat != NULL)
		{
			free(rwCopy.dat);
		}
	}
		
	if (idxU32.dat != NULL)
	{
		free(idxU32.dat);
	}
	if (nck.dat != NULL)
	{
		free(nck.dat);
	}
	if (rs.dat != NULL)
	{
		free(rs.dat);
	}
	if (ncw.dat != NULL)
	{
		free(ncw.dat);
	}
	if (sep1.sep != NULL)
	{
		free(sep1.sep);
	}
	if (sep1.rw != NULL)
	{
		free(sep1.rw);
	}

	return;
}