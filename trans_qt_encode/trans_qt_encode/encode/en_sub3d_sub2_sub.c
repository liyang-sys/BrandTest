#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

void en_sub3d_sub2_sub(Uint8_Dat* z, Uint32_Dat* r, int lB, float pr, Uint32_Dat* nr, float scr, float scrc, Uint32_Dat* cfk, Uint32_Dat* nck, int thd)
{
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
	extern int ptr;
	extern uchar *bin;
	int binBytePos = ptr & 7;

	int typ = 0;
	if (thd % 3 == 0)
	{
		typ = 2;
	}
	else
	{
		typ = 1;
	}
	en_r_cr(r, lB, nr, scr, scrc, pr, z);

	int thdh, thdq;
	if (thd == 2)
	{
		thdh = 1;
		SEP sepCfk = separate0((int *)(cfk->dat), thdh, (int)(cfk->len), 0);
		float cr1 = 0;
		int crc1 = 0;
		Uint32_Dat r1;
		int lB1 = 0;
		float pr = 0;
		Uint32_Dat nr1;
		Uint32_Dat z1;
		z1.dat = sepCfk.sep;
		z1.len = sepCfk.lensepbit;
		k_criterion((void *)&z1, &cr1, &crc1, &r1, &lB1, &pr, &nr1, NULL, 1);

		for (int i = 0; i < nr1.len - 1; i++)
		{
			nr1.dat[i] = nr1.dat[i + 1];
		}
		nr1.len = nr1.len - 1;

		Uint8_Dat z1Copy;//备份
		z1Copy.len = z1.len;
		z1Copy.dat = (unsigned char*)calloc(z1.len / 8 + 1, sizeof(unsigned char));
		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++)
		{
			z1Copy.dat[i] = sepCfk.sep[i];
		}

		//for(int i = 0; i < )
		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++) //将1变为0，将零变为1.
		{
			(sepCfk.sep[i]) ^= 0xff;
		}
		float cr1i = 0;
		int crc1i = 0;
		k_criterion((void*)&z1, &cr1i, &crc1i, NULL, NULL, NULL, NULL, NULL, 1);
		float scr1 = cr1 + cr1i;
		int scrc1 = crc1 + crc1i;
		en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);
		if (sepCfk.sep != NULL)
		{
			free(sepCfk.sep);
		}
		if (sepCfk.rw != NULL)
		{
			free(sepCfk.rw);
		}
		free(nr1.dat);
		free(r1.dat);
		free(z1Copy.dat);
	}
	else if (thd == 3)
	{
		thdh = 2;
		thdq = 1;
		SEP sepCfk = separate0((int *)(cfk->dat), thdh, (int)(cfk->len), 1);
		Uint32_Dat nc11;
		nc11.len = thdh;
		nc11.dat = calloc(nc11.len, sizeof(unsigned int));
		for (int i = 0; i < nc11.len; i++)
		{
			nc11.dat[i] = nck->dat[i];
		}

		float cr1 = 0;
		int crc1 = 0;
		Uint32_Dat r1;
		int lB1 = 0;
		float pr = 0;
		Uint32_Dat nr1;
		Uint32_Dat z1;
		z1.dat = sepCfk.sep;
		z1.len = sepCfk.lensepbit;
		k_criterion((void *)&z1, &cr1, &crc1, &r1, &lB1, &pr, &nr1, NULL, 1);

		for (int i = 0; i < nr1.len - 1; i++)
		{
			nr1.dat[i] = nr1.dat[i + 1];
		}
		//nr1.len = nr1.len - 1;

		//for(int i = 0; i < )
		Uint8_Dat z1Copy;//备份
		z1Copy.len = z1.len;
		z1Copy.dat = (unsigned char*)calloc(z1.len / 8 + 1, sizeof(unsigned char));
		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++)
		{
			z1Copy.dat[i] = sepCfk.sep[i];
		}

		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++) //将1变为0，将零变为1.
		{
			(sepCfk.sep[i]) ^= 0xff;
		}
		float cr1i = 0;
		int crc1i = 0;
		k_criterion((void*)&z1, &cr1i, &crc1i, NULL, NULL, NULL, NULL, NULL, 1);
		float scr1 = cr1 + cr1i;
		int scrc1 = crc1 + crc1i;

		Uint32_Dat cf11;
		cf11.dat = cfk->dat;
		cf11.len = sepCfk.lrk;
		Uint32_Dat cf11Copy;
		cf11Copy.dat = (unsigned int *)calloc(cf11.len, sizeof(unsigned int));
		cf11Copy.len = cf11.len;
		for (int i = 0; i < cf11.len; i++)
		{
			cf11Copy.dat[i] = cf11.dat[i];
		}
		CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);
		if (crTest.sumcr > 12)
		{
			binBytePos = ptr & 7;
			rem.a = 1;//编L2
			rem.a = rem.a << (15 - binBytePos);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr++; binBytePos++;
		}
		else
		{
			ptr++;
		}
		en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);
		if (crTest.sumcr > 12)
		{
			en_r_cr(&(crTest.r), crTest.lB, &(crTest.nr), crTest.sumcr, crTest.sumcrc, crTest.pr, &(crTest.z));
			//en_r_cr(&())
		}
		else
		{
			typ = 1;
			en_2n_code_new(&cf11Copy, &nc11, typ, thdh);

		}

		free(sepCfk.sep);
		free(sepCfk.rw);
		free(crTest.cf1.dat);
		free(crTest.cf2.dat);
		free(crTest.nc1.dat);
		free(crTest.nc2.dat);
		free(crTest.nr.dat);
		free(crTest.r.dat);
		free(cf11Copy.dat);
		free(nc11.dat);
		free(z1Copy.dat);

	}
	else if (thd == 4)
	{
		int L2 = 0;
		thdh = 2;
		thdq = 1;
		SEP sepCfk = separate0((int *)(cfk->dat), thdh, (int)(cfk->len), 1);
		Uint32_Dat cf12;
		cf12.dat = sepCfk.rw;
		cf12.len = sepCfk.lrw;

		Uint32_Dat cf12Copy;
		cf12Copy.dat = (unsigned int *)calloc(cf12.len, sizeof(unsigned int));
		cf12Copy.len = cf12.len;
		for (int i = 0; i < cf12.len; i++)
		{
			cf12Copy.dat[i] = cf12.dat[i];
		}

		Uint32_Dat nc11;
		nc11.len = thdh;
		nc11.dat = calloc(nc11.len, sizeof(unsigned int));
		for (int i = 0; i < nc11.len; i++)
		{
			nc11.dat[i] = nck->dat[i];
		}

		Uint32_Dat nc12;
		nc12.len = nck->len - thdh;
		nc12.dat = calloc(nc12.len, sizeof(unsigned int));
		for (int i = thdh; i < nck->len; i++)
		{
			nc12.dat[i - thdh] = nck->dat[i];
		}


		float cr1 = 0;
		int crc1 = 0;
		Uint32_Dat r1;
		int lB1 = 0;
		float pr = 0;
		Uint32_Dat nr1;
		Uint32_Dat z1;
		z1.dat = sepCfk.sep;
		z1.len = sepCfk.lensepbit;
		k_criterion((void *)&z1, &cr1, &crc1, &r1, &lB1, &pr, &nr1, NULL, 1);

		for (int i = 0; i < nr1.len - 1; i++)
		{
			nr1.dat[i] = nr1.dat[i + 1];
		}
		nr1.len = nr1.len - 1;

		//for(int i = 0; i < )
		Uint8_Dat z1Copy;//备份
		z1Copy.len = z1.len;
		z1Copy.dat = (unsigned char*)calloc(z1.len / 8 + 1, sizeof(unsigned char));
		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++)
		{
			z1Copy.dat[i] = sepCfk.sep[i];
		}

		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++) //将1变为0，将零变为1.
		{
			(sepCfk.sep[i]) ^= 0xff;
		}
		float cr1i = 0;
		int crc1i = 0;
		k_criterion((void*)&z1, &cr1i, &crc1i, NULL, NULL, NULL, NULL, NULL, 1);
		float scr1 = cr1 + cr1i;
		int scrc1 = crc1 + crc1i;

		Uint32_Dat cf11;
		cf11.dat = cfk->dat;
		cf11.len = sepCfk.lrk;
		Uint32_Dat cf11Copy;
		cf11Copy.dat = (unsigned int *)calloc(cf11.len, sizeof(unsigned int));
		cf11Copy.len = cf11.len;
		for (int i = 0; i < cf11.len; i++)
		{
			cf11Copy.dat[i] = cf11.dat[i];
		}
		CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);
		if (crTest.sumcr > 12)
		{
			CR_TEST_RET crTest1 = cr_test(&cf12, thdq, &nc12);
			if (crTest1.sumcr > 12)
			{
				SFC sfc1 = SFcode(3, 3);
				binBytePos = ptr & 7;
				rem.a = sfc1.code;
				rem.a = rem.a << (16 - binBytePos - sfc1.lb);
				bin[ptr >> 3] |= rem.b[1];
				bin[(ptr >> 3) + 1] |= rem.b[0];
				ptr += sfc1.lb;

				en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1
				en_r_cr(&(crTest.r), crTest.lB, &(crTest.nr), crTest.sumcr, crTest.sumcrc, crTest.pr, &(crTest.z));//bin11
				en_r_cr(&(crTest1.r), crTest1.lB, &(crTest1.nr), crTest1.sumcr, crTest1.sumcrc, crTest1.pr, &(crTest1.z));//bin12
			}
			else
			{
				SFC sfc1 = SFcode(2, 3);
				binBytePos = ptr & 7;
				rem.a = sfc1.code;
				rem.a = rem.a << (16 - binBytePos - sfc1.lb);
				bin[ptr >> 3] |= rem.b[1];
				bin[(ptr >> 3) + 1] |= rem.b[0];
				ptr += sfc1.lb;

				en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1
				en_r_cr(&(crTest.r), crTest.lB, &(crTest.nr), crTest.sumcr, crTest.sumcrc, crTest.pr, &(crTest.z));//bin11
				en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin12
			}
			//free(crTest1.cf1.dat);
			free(crTest1.cf2.dat);
			free(crTest1.nc1.dat);
			free(crTest1.nc2.dat);
			free(crTest1.nr.dat);
			free(crTest1.r.dat);
		}
		else
		{
			SFC sfc1 = SFcode(1, 3);
			binBytePos = ptr & 7;
			rem.a = sfc1.code;
			rem.a = rem.a << (16 - binBytePos - sfc1.lb);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr += sfc1.lb;

			en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1
			en_2n_code_new(&cf11Copy, &nc11, typ, thdh);//bin11
			en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin12
		}

		free(sepCfk.sep);
		free(sepCfk.rw);
		free(crTest.cf1.dat);
		free(crTest.cf2.dat);
		free(crTest.nc1.dat);
		free(crTest.nc2.dat);
		free(crTest.nr.dat);
		free(crTest.r.dat);
		free(cf11Copy.dat);
		free(cf12Copy.dat);
		free(nc11.dat);
		free(nc12.dat);
		free(z1Copy.dat);
	}
	else if (thd == 6)
	{
		thdh = 3;
		thdq = 2;
		int thdt = 1, typ2 = 1, L2 = 0, L3 = 0;
		SEP sepCfk = separate0((int *)(cfk->dat), thdh, (int)(cfk->len), 1);
		Uint32_Dat cf12;
		cf12.dat = sepCfk.rw;
		cf12.len = sepCfk.lrw;

		Uint32_Dat cf12Copy;
		cf12Copy.dat = (unsigned int *)calloc(cf12.len, sizeof(unsigned int));
		cf12Copy.len = cf12.len;
		for (int i = 0; i < cf12.len; i++)
		{
			cf12Copy.dat[i] = cf12.dat[i];
		}

		Uint32_Dat nc11;
		nc11.len = thdh;
		nc11.dat = calloc(nc11.len, sizeof(unsigned int));
		for (int i = 0; i < nc11.len; i++)
		{
			nc11.dat[i] = nck->dat[i];
		}

		Uint32_Dat nc12;
		nc12.len = nck->len - thdh;
		nc12.dat = calloc(nc12.len, sizeof(unsigned int));
		for (int i = thdh; i < nck->len; i++)
		{
			nc12.dat[i - thdh] = nck->dat[i];
		}


		float cr1 = 0;
		int crc1 = 0;
		Uint32_Dat r1;
		int lB1 = 0;
		float pr = 0;
		Uint32_Dat nr1;
		Uint32_Dat z1;
		z1.dat = sepCfk.sep;
		z1.len = sepCfk.lensepbit;
		k_criterion((void *)&z1, &cr1, &crc1, &r1, &lB1, &pr, &nr1, NULL, 1);

		for (int i = 0; i < nr1.len - 1; i++)
		{
			nr1.dat[i] = nr1.dat[i + 1];
		}
		nr1.len = nr1.len - 1;

		//for(int i = 0; i < )
		Uint8_Dat z1Copy;//备份
		z1Copy.len = z1.len;
		z1Copy.dat = (unsigned char*)calloc(z1.len / 8 + 1, sizeof(unsigned char));
		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++)
		{
			z1Copy.dat[i] = sepCfk.sep[i];
		}

		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++) //将1变为0，将零变为1.
		{
			(sepCfk.sep[i]) ^= 0xff;
		}
		float cr1i = 0;
		int crc1i = 0;
		k_criterion((void*)&z1, &cr1i, &crc1i, NULL, NULL, NULL, NULL, NULL, 1);
		float scr1 = cr1 + cr1i;
		int scrc1 = crc1 + crc1i;

		Uint32_Dat cf11;
		cf11.dat = cfk->dat;
		cf11.len = sepCfk.lrk;
		Uint32_Dat cf11Copy;
		cf11Copy.dat = (unsigned int *)calloc(cf11.len, sizeof(unsigned int));
		cf11Copy.len = cf11.len;
		for (int i = 0; i < cf11.len; i++)
		{
			cf11Copy.dat[i] = cf11.dat[i];
		}


		CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);
		if (crTest.sumcr > 12)
		{
			L2 = 1;
			///////
			Uint32_Dat crTestCf1Copy;//调用之后数据会被改变，所以需要拷贝
			crTestCf1Copy.len = crTest.cf1.len;
			crTestCf1Copy.dat = calloc(crTest.cf1.len, sizeof(unsigned int));
			for (int i = 0; i < crTest.cf1.len; i++)
			{
				crTestCf1Copy.dat[i] = crTest.cf1.dat[i];
			}
			CR_TEST_RET crTest1 = cr_test(&crTestCf1Copy, thdt, &crTest.nc1);
			//free(crTestCf1Copy.dat);

			//CR_TEST_RET crTest1 = cr_test(&crTest.cf1, thdt, &crTest.nc1);
			if (crTest1.sumcr > 12)
			{
				L3 = 1;
			}
			CR_TEST_RET crTest2 = cr_test(&cf12, thdq, &nc12);
			if (crTest2.sumcr > 12)
			{
				L2 = 2;

				Uint32_Dat crTest2Cf1Copy;//调用之后数据会被改变，所以需要拷贝
				crTest2Cf1Copy.len = crTest2.cf1.len;
				crTest2Cf1Copy.dat = calloc(crTest2.cf1.len, sizeof(unsigned int));
				for (int i = 0; i < crTest2.cf1.len; i++)
				{
					crTest2Cf1Copy.dat[i] = crTest2.cf1.dat[i];
				}
				CR_TEST_RET crTest3 = cr_test(&crTest2Cf1Copy, thdt, &crTest2.nc1);
				//free(crTest2Cf1Copy.dat);

				//CR_TEST_RET crTest3 = cr_test(&crTest2.cf1, thdt, &crTest2.nc1);
				if ((crTest3.sumcr > 16) && (L3 == 1))
				{
					L3 = 2;
					SFC sfc1 = SFcode(L2 + 1, 3);
					binBytePos = ptr & 7;
					rem.a = sfc1.code;
					rem.a = rem.a << (16 - binBytePos - sfc1.lb);
					bin[ptr >> 3] |= rem.b[1];
					bin[(ptr >> 3) + 1] |= rem.b[0];
					ptr += sfc1.lb;

					SFC sfc2 = SFcode(L3 + 1, 3);
					binBytePos = ptr & 7;
					rem.a = sfc2.code;
					rem.a = rem.a << (16 - binBytePos - sfc2.lb);
					bin[ptr >> 3] |= rem.b[1];
					bin[(ptr >> 3) + 1] |= rem.b[0];
					ptr += sfc2.lb;

					en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1
					en_r_cr(&(crTest.r), crTest.lB, &(crTest.nr), crTest.sumcr, crTest.sumcrc, crTest.pr, &(crTest.z));//bin_r2A
					if (crTest1.sumcr > 12)
					{
						en_r_cr(&(crTest1.r), crTest1.lB, &(crTest1.nr), crTest1.sumcr, crTest1.sumcrc, crTest1.pr, &(crTest1.z));//bin21
					}
					else
					{
						en_2n_code_new(&crTest.cf1, &crTest.nc1, typ2, thdq);//bin21
					}
					en_r_cr(&(crTest2.r), crTest2.lB, &(crTest2.nr), crTest2.sumcr, crTest2.sumcrc, crTest2.pr, &(crTest2.z));//bin_r2B
					en_r_cr(&(crTest3.r), crTest3.lB, &(crTest3.nr), crTest3.sumcr, crTest3.sumcrc, crTest3.pr, &(crTest3.z));//bin23

				}
				else
				{
					SFC sfc1 = SFcode(L2 + 1, 3);
					binBytePos = ptr & 7;
					rem.a = sfc1.code;
					rem.a = rem.a << (16 - binBytePos - sfc1.lb);
					bin[ptr >> 3] |= rem.b[1];
					bin[(ptr >> 3) + 1] |= rem.b[0];
					ptr += sfc1.lb;

					SFC sfc2 = SFcode(L3 + 1, 3);
					binBytePos = ptr & 7;
					rem.a = sfc2.code;
					rem.a = rem.a << (16 - binBytePos - sfc2.lb);
					bin[ptr >> 3] |= rem.b[1];
					bin[(ptr >> 3) + 1] |= rem.b[0];
					ptr += sfc2.lb;

					en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1
					en_r_cr(&(crTest.r), crTest.lB, &(crTest.nr), crTest.sumcr, crTest.sumcrc, crTest.pr, &(crTest.z));//bin_r2A
					if (crTest1.sumcr > 12)
					{
						en_r_cr(&(crTest1.r), crTest1.lB, &(crTest1.nr), crTest1.sumcr, crTest1.sumcrc, crTest1.pr, &(crTest1.z));//bin21
					}
					else
					{
						en_2n_code_new(&crTest.cf1, &crTest.nc1, typ2, thdq);//bin21
					}
					en_r_cr(&(crTest2.r), crTest2.lB, &(crTest2.nr), crTest2.sumcr, crTest2.sumcrc, crTest2.pr, &(crTest2.z));//bin_r2B
					en_2n_code_new(&crTest2.cf1, &crTest2.nc1, typ2, thdq);//bin23
				}
				//free(crTest3.cf1.dat);
				free(crTest3.cf2.dat);
				free(crTest3.nc1.dat);
				free(crTest3.nc2.dat);
				free(crTest3.nr.dat);
				free(crTest3.r.dat);
			}
			else
			{
				SFC sfc1 = SFcode(L2 + 1, 3);
				binBytePos = ptr & 7;
				rem.a = sfc1.code;
				rem.a = rem.a << (16 - binBytePos - sfc1.lb);
				bin[ptr >> 3] |= rem.b[1];
				bin[(ptr >> 3) + 1] |= rem.b[0];
				ptr += sfc1.lb;

				SFC sfc2 = SFcode(L3 + 1, 3);
				binBytePos = ptr & 7;
				rem.a = sfc2.code;
				rem.a = rem.a << (16 - binBytePos - sfc2.lb);
				bin[ptr >> 3] |= rem.b[1];
				bin[(ptr >> 3) + 1] |= rem.b[0];
				ptr += sfc2.lb;

				en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1
				en_r_cr(&(crTest.r), crTest.lB, &(crTest.nr), crTest.sumcr, crTest.sumcrc, crTest.pr, &(crTest.z));//bin_r2A
				if (crTest1.sumcr > 12)
				{
					en_r_cr(&(crTest1.r), crTest1.lB, &(crTest1.nr), crTest1.sumcr, crTest1.sumcrc, crTest1.pr, &(crTest1.z));//bin21
				}
				else
				{
					en_2n_code_new(&crTest.cf1, &crTest.nc1, typ2, thdq);//bin21
				}

				//FILE* fp = fopen("cf12Copy.txt", "wb");
				//fwrite(cf12Copy.dat, sizeof(unsigned int), cf12Copy.len, fp);
				//fclose(fp);

				en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin12
			}

			//free(crTest1.cf1.dat);
			free(crTest1.cf2.dat);
			free(crTest1.nc1.dat);
			free(crTest1.nc2.dat);
			free(crTest1.nr.dat);
			free(crTest1.r.dat);

			//free(crTest2.cf1.dat);
			free(crTest2.cf2.dat);
			free(crTest2.nc1.dat);
			free(crTest2.nc2.dat);
			free(crTest2.nr.dat);
			free(crTest2.r.dat);
		}
		else
		{
			SFC sfc1 = SFcode(L2 + 1, 3);
			binBytePos = ptr & 7;
			rem.a = sfc1.code;
			rem.a = rem.a << (16 - binBytePos - sfc1.lb);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr += sfc1.lb;

			SFC sfc2 = SFcode(L3 + 1, 3);
			binBytePos = ptr & 7;
			rem.a = sfc2.code;
			rem.a = rem.a << (16 - binBytePos - sfc2.lb);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr += sfc2.lb;

			en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1
			en_2n_code_new(&cf11Copy, &nc11, typ, thdh);//bin11
			en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin12
		}

		free(sepCfk.sep);
		free(sepCfk.rw);
		free(crTest.cf1.dat);
		free(crTest.cf2.dat);
		free(crTest.nc1.dat);
		free(crTest.nc2.dat);
		free(crTest.nr.dat);
		free(crTest.r.dat);
		free(cf11Copy.dat);
		free(cf12Copy.dat);
		free(nc11.dat);
		free(nc12.dat);
		free(z1Copy.dat);
	}
	else if (thd == 8)
	{
		thdh = 4;
		thdq = 2;
		int thdt = 1, typ2 = 1, L2 = 0, L3 = 0;
		SEP sepCfk = separate0((int *)(cfk->dat), thdh, (int)(cfk->len), 1);
		Uint32_Dat cf12;
		cf12.dat = sepCfk.rw;
		cf12.len = sepCfk.lrw;

		Uint32_Dat cf12Copy;
		cf12Copy.dat = (unsigned int *)calloc(cf12.len, sizeof(unsigned int));
		cf12Copy.len = cf12.len;
		for (int i = 0; i < cf12.len; i++)
		{
			cf12Copy.dat[i] = cf12.dat[i];
		}

		Uint32_Dat nc11;
		nc11.len = thdh;
		nc11.dat = calloc(nc11.len, sizeof(unsigned int));
		for (int i = 0; i < nc11.len; i++)
		{
			nc11.dat[i] = nck->dat[i];
		}

		Uint32_Dat nc12;
		nc12.len = nck->len - thdh;
		nc12.dat = calloc(nc12.len, sizeof(unsigned int));
		for (int i = thdh; i < nck->len; i++)
		{
			nc12.dat[i - thdh] = nck->dat[i];
		}


		float cr1 = 0;
		int crc1 = 0;
		Uint32_Dat r1;
		int lB1 = 0;
		float pr = 0;
		Uint32_Dat nr1;
		Uint32_Dat z1;
		z1.dat = sepCfk.sep;
		z1.len = sepCfk.lensepbit;
		k_criterion((void *)&z1, &cr1, &crc1, &r1, &lB1, &pr, &nr1, NULL, 1);

		for (int i = 0; i < nr1.len - 1; i++)
		{
			nr1.dat[i] = nr1.dat[i + 1];
		}
		nr1.len = nr1.len - 1;

		//for(int i = 0; i < )
		Uint8_Dat z1Copy;//备份
		z1Copy.len = z1.len;
		z1Copy.dat = (unsigned char*)calloc(z1.len / 8 + 1, sizeof(unsigned char));
		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++)
		{
			z1Copy.dat[i] = sepCfk.sep[i];
		}

		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++) //将1变为0，将零变为1.
		{
			(sepCfk.sep[i]) ^= 0xff;
		}
		float cr1i = 0;
		int crc1i = 0;
		k_criterion((void*)&z1, &cr1i, &crc1i, NULL, NULL, NULL, NULL, NULL, 1);
		float scr1 = cr1 + cr1i;
		int scrc1 = crc1 + crc1i;

		Uint32_Dat cf11;
		cf11.dat = cfk->dat;
		cf11.len = sepCfk.lrk;
		Uint32_Dat cf11Copy;
		cf11Copy.dat = (unsigned int *)calloc(cf11.len, sizeof(unsigned int));
		cf11Copy.len = cf11.len;
		for (int i = 0; i < cf11.len; i++)
		{
			cf11Copy.dat[i] = cf11.dat[i];
		}


		//CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);
		//CR_TEST_RET crTest1 = cr_test(&crTest.cf1, thdt, &crTest.nc1);
		//CR_TEST_RET crTest2 = cr_test(&crTest.cf2, thdt, &crTest.nc2);
		//CR_TEST_RET crTest3 = cr_test(&cf12, thdq, &nc12);
		//CR_TEST_RET crTest4 = cr_test(&crTest3.cf1, thdt, &crTest3.nc1);
		//CR_TEST_RET crTest5 = cr_test(&crTest3.cf2, thdt, &crTest3.nc2);
		CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);

		Uint32_Dat crTestCf1Copy;//调用之后数据会被改变，所以需要拷贝
		crTestCf1Copy.len = crTest.cf1.len;
		crTestCf1Copy.dat = calloc(crTest.cf1.len, sizeof(unsigned int));
		for (int i = 0; i < crTest.cf1.len; i++)
		{
			crTestCf1Copy.dat[i] = crTest.cf1.dat[i];
		}
		CR_TEST_RET crTest1 = cr_test(&crTestCf1Copy, thdt, &crTest.nc1);
		//free(crTestCf1Copy.dat);

		Uint32_Dat crTestCf2Copy;//调用之后数据会被改变，所以需要拷贝
		crTestCf2Copy.len = crTest.cf2.len;
		crTestCf2Copy.dat = calloc(crTest.cf2.len, sizeof(unsigned int));
		for (int i = 0; i < crTest.cf2.len; i++)
		{
			crTestCf2Copy.dat[i] = crTest.cf2.dat[i];
		}
		CR_TEST_RET crTest2 = cr_test(&crTestCf2Copy, thdt, &crTest.nc2);
		//free(crTestCf2Copy.dat);

		CR_TEST_RET crTest3 = cr_test(&cf12, thdq, &nc12);

		Uint32_Dat crTest3Cf1Copy;//调用之后数据会被改变，所以需要拷贝
		crTest3Cf1Copy.len = crTest3.cf1.len;
		crTest3Cf1Copy.dat = calloc(crTest3.cf1.len, sizeof(unsigned int));
		for (int i = 0; i < crTest3.cf1.len; i++)
		{
			crTest3Cf1Copy.dat[i] = crTest3.cf1.dat[i];
		}
		CR_TEST_RET crTest4 = cr_test(&crTest3Cf1Copy, thdt, &crTest3.nc1);
		//free(crTest3Cf1Copy.dat);

		Uint32_Dat crTest3Cf2Copy;//调用之后数据会被改变，所以需要拷贝
		crTest3Cf2Copy.len = crTest3.cf2.len;
		crTest3Cf2Copy.dat = calloc(crTest3.cf2.len, sizeof(unsigned int));
		for (int i = 0; i < crTest3.cf2.len; i++)
		{
			crTest3Cf2Copy.dat[i] = crTest3.cf2.dat[i];
		}
		CR_TEST_RET crTest5 = cr_test(&crTest3Cf2Copy, thdt, &crTest3.nc2);
		//free(crTest3Cf2Copy.dat);
		if (crTest.sumcr > 12)
		{
			L2 = 1;

			if (crTest1.sumcr > 12)
			{
				L3 = 1;
				if (crTest2.sumcr > 12)
				{
					L3 = 2;
				}

			}
			if (crTest3.sumcr > 12)
			{
				L2 = 2;

				if ((crTest4.sumcr > 16) && (L3 == 2))
				{
					L3 = 3;

					if (crTest5.sumcr > 16)
					{
						L3 = 4;
					}

				}
			}
		}

		SFC sfc1 = SFcode(L2 + 1, 3);
		binBytePos = ptr & 7;
		rem.a = sfc1.code;
		rem.a = rem.a << (16 - binBytePos - sfc1.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc1.lb;

		SFC sfc2 = SFcode(L3 + 1, 5);
		binBytePos = ptr & 7;
		rem.a = sfc2.code;
		rem.a = rem.a << (16 - binBytePos - sfc2.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc2.lb;

		en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1

		if (crTest.sumcr > 12)
		{
			//L2 = 1;
			en_r_cr(&(crTest.r), crTest.lB, &(crTest.nr), crTest.sumcr, crTest.sumcrc, crTest.pr, &(crTest.z));//bin_r2A
			if (crTest1.sumcr > 12)
			{
				//L3 = 1;
				en_r_cr(&(crTest1.r), crTest1.lB, &(crTest1.nr), crTest1.sumcr, crTest1.sumcrc, crTest1.pr, &(crTest1.z));//bin21
				if (crTest2.sumcr > 12)
				{
					//L3 = 2;
					en_r_cr(&(crTest2.r), crTest2.lB, &(crTest2.nr), crTest2.sumcr, crTest2.sumcrc, crTest2.pr, &(crTest2.z));//bin22
				}
				else
				{
					en_2n_code_new(&crTest.cf2, &crTest.nc2, typ, thdq);//bin22
				}
			}
			else
			{
				en_2n_code_new(&crTest.cf1, &crTest.nc1, typ, thdq);//bin21
				en_2n_code_new(&crTest.cf2, &crTest.nc2, typ, thdq);//bin22
			}

			if (crTest3.sumcr > 12)
			{
				//L2 = 2;
				en_r_cr(&(crTest3.r), crTest3.lB, &(crTest3.nr), crTest3.sumcr, crTest3.sumcrc, crTest3.pr, &(crTest3.z));//bin_r2B

				if ((crTest4.sumcr > 16) && (L3 == 2))
				{
					//L3 = 3;
					en_r_cr(&(crTest4.r), crTest4.lB, &(crTest4.nr), crTest4.sumcr, crTest4.sumcrc, crTest4.pr, &(crTest4.z));//bin23
					if (crTest5.sumcr > 16)
					{
						//L3 = 4;
						en_r_cr(&(crTest5.r), crTest5.lB, &(crTest5.nr), crTest5.sumcr, crTest5.sumcrc, crTest5.pr, &(crTest5.z));//bin24
					}
					else
					{
						en_2n_code_new(&crTest3.cf2, &crTest3.nc2, typ, thdq);//bin24
					}

				}
				else
				{
					en_2n_code_new(&crTest3.cf1, &crTest3.nc1, typ, thdq);//bin21
					en_2n_code_new(&crTest3.cf2, &crTest3.nc2, typ, thdq);//bin22
				}

			}
			else
			{
				en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin12
			}


		}
		else
		{
			en_2n_code_new(&cf11Copy, &nc11, typ, thdh);//bin23
			en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin24
		}


		//CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);
		//CR_TEST_RET crTest1 = cr_test(&crTest.cf1, thdt, &crTest.nc1);
		//CR_TEST_RET crTest2 = cr_test(&crTest.cf2, thdt, &crTest.nc2);
		//CR_TEST_RET crTest3 = cr_test(&cf12, thdq, &nc12);
		//CR_TEST_RET crTest4 = cr_test(&crTest3.cf1, thdt, &crTest3.nc1);
		//CR_TEST_RET crTest5 = cr_test(&crTest3.cf2, thdt, &crTest3.nc2);


		free(sepCfk.sep);
		free(sepCfk.rw);
		free(crTest.cf1.dat);
		free(crTest.cf2.dat);
		free(crTest.nc1.dat);
		free(crTest.nc2.dat);
		free(crTest.nr.dat);
		free(crTest.r.dat);
		free(cf11Copy.dat);
		free(cf12Copy.dat);
		free(nc11.dat);
		free(nc12.dat);
		free(z1Copy.dat);


		free(crTest1.cf1.dat);
		free(crTest1.cf2.dat);
		free(crTest1.nc1.dat);
		free(crTest1.nc2.dat);
		free(crTest1.nr.dat);
		free(crTest1.r.dat);

		free(crTest2.cf1.dat);
		free(crTest2.cf2.dat);
		free(crTest2.nc1.dat);
		free(crTest2.nc2.dat);
		free(crTest2.nr.dat);
		free(crTest2.r.dat);

		//free(crTest3.cf1.dat);
		free(crTest3.cf2.dat);
		free(crTest3.nc1.dat);
		free(crTest3.nc2.dat);
		free(crTest3.nr.dat);
		free(crTest3.r.dat);

		free(crTest4.cf1.dat);
		free(crTest4.cf2.dat);
		free(crTest4.nc1.dat);
		free(crTest4.nc2.dat);
		free(crTest4.nr.dat);
		free(crTest4.r.dat);

		free(crTest5.cf1.dat);
		free(crTest5.cf2.dat);
		free(crTest5.nc1.dat);
		free(crTest5.nc2.dat);
		free(crTest5.nr.dat);
		free(crTest5.r.dat);
	}
	else if (thd == 12)
	{
		thdh = 6;
		thdq = 3;
		int thdt = 2, typ3 = 1, L2 = 0, L3 = 0;
		SEP sepCfk = separate0((int *)(cfk->dat), thdh, (int)(cfk->len), 1);
		Uint32_Dat cf12;
		cf12.dat = sepCfk.rw;
		cf12.len = sepCfk.lrw;

		Uint32_Dat cf12Copy;
		cf12Copy.dat = (unsigned int *)calloc(cf12.len, sizeof(unsigned int));
		cf12Copy.len = cf12.len;
		for (int i = 0; i < cf12.len; i++)
		{
			cf12Copy.dat[i] = cf12.dat[i];
		}

		Uint32_Dat nc11;
		nc11.len = thdh;
		nc11.dat = calloc(nc11.len, sizeof(unsigned int));
		for (int i = 0; i < nc11.len; i++)
		{
			nc11.dat[i] = nck->dat[i];
		}

		Uint32_Dat nc12;
		nc12.len = nck->len - thdh;
		nc12.dat = calloc(nc12.len, sizeof(unsigned int));
		for (int i = thdh; i < nck->len; i++)
		{
			nc12.dat[i - thdh] = nck->dat[i];
		}


		float cr1 = 0;
		int crc1 = 0;
		Uint32_Dat r1;
		int lB1 = 0;
		float pr = 0;
		Uint32_Dat nr1;
		Uint32_Dat z1;
		z1.dat = sepCfk.sep;
		z1.len = sepCfk.lensepbit;
		k_criterion((void *)&z1, &cr1, &crc1, &r1, &lB1, &pr, &nr1, NULL, 1);

		for (int i = 0; i < nr1.len - 1; i++)
		{
			nr1.dat[i] = nr1.dat[i + 1];
		}
		nr1.len = nr1.len - 1;

		//for(int i = 0; i < )
		Uint8_Dat z1Copy;//备份
		z1Copy.len = z1.len;
		z1Copy.dat = (unsigned char*)calloc(z1.len / 8 + 1, sizeof(unsigned char));
		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++)
		{
			z1Copy.dat[i] = sepCfk.sep[i];
		}

		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++) //将1变为0，将零变为1.
		{
			(sepCfk.sep[i]) ^= 0xff;
		}
		float cr1i = 0;
		int crc1i = 0;
		k_criterion((void*)&z1, &cr1i, &crc1i, NULL, NULL, NULL, NULL, NULL, 1);
		float scr1 = cr1 + cr1i;
		int scrc1 = crc1 + crc1i;

		Uint32_Dat cf11;
		cf11.dat = cfk->dat;
		cf11.len = sepCfk.lrk;
		Uint32_Dat cf11Copy;
		cf11Copy.dat = (unsigned int *)calloc(cf11.len, sizeof(unsigned int));
		cf11Copy.len = cf11.len;
		for (int i = 0; i < cf11.len; i++)
		{
			cf11Copy.dat[i] = cf11.dat[i];
		}


		CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);

		Uint32_Dat crTestCf1Copy;//调用之后数据会被改变，所以需要拷贝
		crTestCf1Copy.len = crTest.cf1.len;
		crTestCf1Copy.dat = calloc(crTest.cf1.len, sizeof(unsigned int));
		for (int i = 0; i < crTest.cf1.len; i++)
		{
			crTestCf1Copy.dat[i] = crTest.cf1.dat[i];
		}
		CR_TEST_RET crTest1 = cr_test(&crTestCf1Copy, thdt, &crTest.nc1);
		//free(crTestCf1Copy.dat);

		Uint32_Dat crTestCf2Copy;//调用之后数据会被改变，所以需要拷贝
		crTestCf2Copy.len = crTest.cf2.len;
		crTestCf2Copy.dat = calloc(crTest.cf2.len, sizeof(unsigned int));
		for (int i = 0; i < crTest.cf2.len; i++)
		{
			crTestCf2Copy.dat[i] = crTest.cf2.dat[i];
		}
		CR_TEST_RET crTest2 = cr_test(&crTestCf2Copy, thdt, &crTest.nc2);
		//free(crTestCf2Copy.dat);

		CR_TEST_RET crTest3 = cr_test(&cf12, thdq, &nc12);

		Uint32_Dat crTest3Cf1Copy;//调用之后数据会被改变，所以需要拷贝
		crTest3Cf1Copy.len = crTest3.cf1.len;
		crTest3Cf1Copy.dat = calloc(crTest3.cf1.len, sizeof(unsigned int));
		for (int i = 0; i < crTest3.cf1.len; i++)
		{
			crTest3Cf1Copy.dat[i] = crTest3.cf1.dat[i];
		}
		CR_TEST_RET crTest4 = cr_test(&crTest3Cf1Copy, thdt, &crTest3.nc1);
		//free(crTest3Cf1Copy.dat);

		Uint32_Dat crTest3Cf2Copy;//调用之后数据会被改变，所以需要拷贝
		crTest3Cf2Copy.len = crTest3.cf2.len;
		crTest3Cf2Copy.dat = calloc(crTest3.cf2.len, sizeof(unsigned int));
		for (int i = 0; i < crTest3.cf2.len; i++)
		{
			crTest3Cf2Copy.dat[i] = crTest3.cf2.dat[i];
		}
		CR_TEST_RET crTest5 = cr_test(&crTest3Cf2Copy, thdt, &crTest3.nc2);
		//free(crTest3Cf2Copy.dat);

		if (crTest.sumcr > 12)
		{
			L2 = 1;

			if (crTest1.sumcr > 12)
			{
				L3 = 1;
				if (crTest2.sumcr > 12)
				{
					L3 = 2;
				}

			}
			if (crTest3.sumcr > 12)
			{
				L2 = 2;

				if ((crTest4.sumcr > 16) && (L3 == 2))
				{
					L3 = 3;

					if (crTest5.sumcr > 16)
					{
						L3 = 4;
					}

				}
			}
		}

		SFC sfc1 = SFcode(L2 + 1, 3);
		binBytePos = ptr & 7;
		rem.a = sfc1.code;
		rem.a = rem.a << (16 - binBytePos - sfc1.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc1.lb;

		SFC sfc2 = SFcode(L3 + 1, 5);
		binBytePos = ptr & 7;
		rem.a = sfc2.code;
		rem.a = rem.a << (16 - binBytePos - sfc2.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc2.lb;

		en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1

		if (crTest.sumcr > 12)
		{
			//L2 = 1;
			en_r_cr(&(crTest.r), crTest.lB, &(crTest.nr), crTest.sumcr, crTest.sumcrc, crTest.pr, &(crTest.z));//bin_r2A
			if (crTest1.sumcr > 12)
			{
				//L3 = 1;
				en_r_cr(&(crTest1.r), crTest1.lB, &(crTest1.nr), crTest1.sumcr, crTest1.sumcrc, crTest1.pr, &(crTest1.z));//bin_r3A
				en_2n_code_new(&crTest1.cf1, &crTest1.nc1, typ3, thdt);//bin31
				if (crTest2.sumcr > 12)
				{
					//L3 = 2;
					en_r_cr(&(crTest2.r), crTest2.lB, &(crTest2.nr), crTest2.sumcr, crTest2.sumcrc, crTest2.pr, &(crTest2.z));//bin_r3B
					en_2n_code_new(&crTest2.cf1, &crTest2.nc1, typ3, thdt);//bin33
				}
				else
				{
					en_2n_code_new(&crTest.cf2, &crTest.nc2, typ, thdq);//bin22
				}
			}
			else
			{
				en_2n_code_new(&crTest.cf1, &crTest.nc1, typ, thdq);//bin21
				en_2n_code_new(&crTest.cf2, &crTest.nc2, typ, thdq);//bin22
			}

			if (crTest3.sumcr > 12)
			{
				//L2 = 2;
				en_r_cr(&(crTest3.r), crTest3.lB, &(crTest3.nr), crTest3.sumcr, crTest3.sumcrc, crTest3.pr, &(crTest3.z));//bin_r2B

				if ((crTest4.sumcr > 16) && (L3 == 2))
				{
					//L3 = 3;
					en_r_cr(&(crTest4.r), crTest4.lB, &(crTest4.nr), crTest4.sumcr, crTest4.sumcrc, crTest4.pr, &(crTest4.z));//bin_r3C
					en_2n_code_new(&crTest4.cf1, &crTest4.nc1, typ3, thdt);//bin35
					if (crTest5.sumcr > 16)
					{
						//L3 = 4;
						en_r_cr(&(crTest5.r), crTest5.lB, &(crTest5.nr), crTest5.sumcr, crTest5.sumcrc, crTest5.pr, &(crTest5.z));//bin24
						en_2n_code_new(&crTest5.cf1, &crTest5.nc1, typ3, thdt);//bin37
					}
					else
					{
						en_2n_code_new(&crTest3.cf2, &crTest3.nc2, typ, thdq);//bin24
					}

				}
				else
				{
					en_2n_code_new(&crTest3.cf1, &crTest3.nc1, typ, thdq);//bin21
					en_2n_code_new(&crTest3.cf2, &crTest3.nc2, typ, thdq);//bin22
				}

			}
			else
			{
				en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin12
			}


		}
		else
		{
			en_2n_code_new(&cf11Copy, &nc11, typ, thdh);//bin11
			en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin12
		}


		//CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);
		//CR_TEST_RET crTest1 = cr_test(&crTest.cf1, thdt, &crTest.nc1);
		//CR_TEST_RET crTest2 = cr_test(&crTest.cf2, thdt, &crTest.nc2);
		//CR_TEST_RET crTest3 = cr_test(&cf12, thdq, &nc12);
		//CR_TEST_RET crTest4 = cr_test(&crTest3.cf1, thdt, &crTest3.nc1);
		//CR_TEST_RET crTest5 = cr_test(&crTest3.cf2, thdt, &crTest3.nc2);


		free(sepCfk.sep);
		free(sepCfk.rw);
		free(crTest.cf1.dat);
		free(crTest.cf2.dat);
		free(crTest.nc1.dat);
		free(crTest.nc2.dat);
		free(crTest.nr.dat);
		free(crTest.r.dat);
		free(cf11Copy.dat);
		free(cf12Copy.dat);
		free(nc11.dat);
		free(nc12.dat);
		free(z1Copy.dat);


		free(crTest1.cf1.dat);
		free(crTest1.cf2.dat);
		free(crTest1.nc1.dat);
		free(crTest1.nc2.dat);
		free(crTest1.nr.dat);
		free(crTest1.r.dat);

		free(crTest2.cf1.dat);
		free(crTest2.cf2.dat);
		free(crTest2.nc1.dat);
		free(crTest2.nc2.dat);
		free(crTest2.nr.dat);
		free(crTest2.r.dat);

		//free(crTest3.cf1.dat);
		free(crTest3.cf2.dat);
		free(crTest3.nc1.dat);
		free(crTest3.nc2.dat);
		free(crTest3.nr.dat);
		free(crTest3.r.dat);

		free(crTest4.cf1.dat);
		free(crTest4.cf2.dat);
		free(crTest4.nc1.dat);
		free(crTest4.nc2.dat);
		free(crTest4.nr.dat);
		free(crTest4.r.dat);

		free(crTest5.cf1.dat);
		free(crTest5.cf2.dat);
		free(crTest5.nc1.dat);
		free(crTest5.nc2.dat);
		free(crTest5.nr.dat);
		free(crTest5.r.dat);
	}
	else if (thd >= 16)
	{
	    thdh = thd / 2;
		thdq = thd / 4; 
		int thdt = thd / 8;
		int  L2 = 0, L3 = 0;
		SEP sepCfk = separate0((int *)(cfk->dat), thdh, (int)(cfk->len), 1);
		Uint32_Dat cf12;
		cf12.dat = sepCfk.rw;
		cf12.len = sepCfk.lrw;

		Uint32_Dat cf12Copy;
		cf12Copy.dat = (unsigned int *)calloc(cf12.len, sizeof(unsigned int));
		cf12Copy.len = cf12.len;
		for (int i = 0; i < cf12.len; i++)
		{
			cf12Copy.dat[i] = cf12.dat[i];
		}

		Uint32_Dat nc11;
		nc11.len = thdh;
		nc11.dat = calloc(nc11.len, sizeof(unsigned int));
		for (int i = 0; i < nc11.len; i++)
		{
			nc11.dat[i] = nck->dat[i];
		}

		Uint32_Dat nc12;
		nc12.len = nck->len - thdh;
		nc12.dat = calloc(nc12.len, sizeof(unsigned int));
		for (int i = thdh; i < nck->len; i++)
		{
			nc12.dat[i - thdh] = nck->dat[i];
		}


		float cr1 = 0;
		int crc1 = 0;
		Uint32_Dat r1;
		int lB1 = 0;
		float pr = 0;
		Uint32_Dat nr1;
		Uint32_Dat z1;
		z1.dat = sepCfk.sep;
		z1.len = sepCfk.lensepbit;
		k_criterion((void *)&z1, &cr1, &crc1, &r1, &lB1, &pr, &nr1, NULL, 1);

		for (int i = 0; i < nr1.len - 1; i++)
		{
			nr1.dat[i] = nr1.dat[i + 1];
		}
		nr1.len = nr1.len - 1;

		//for(int i = 0; i < )
		Uint8_Dat z1Copy;//备份
		z1Copy.len = z1.len;
		z1Copy.dat = (unsigned char*)calloc(z1.len / 8 + 1, sizeof(unsigned char));
		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++)
		{
			z1Copy.dat[i] = sepCfk.sep[i];
		}

		for (int i = 0; i < (((sepCfk.lensepbit) >> 3) + 1); i++) //将1变为0，将零变为1.
		{
			(sepCfk.sep[i]) ^= 0xff;
		}
		float cr1i = 0;
		int crc1i = 0;
		k_criterion((void*)&z1, &cr1i, &crc1i, NULL, NULL, NULL, NULL, NULL, 1);
		float scr1 = cr1 + cr1i;
		int scrc1 = crc1 + crc1i;

		Uint32_Dat cf11;
		cf11.dat = cfk->dat;
		cf11.len = sepCfk.lrk;
		Uint32_Dat cf11Copy;
		cf11Copy.dat = (unsigned int *)calloc(cf11.len, sizeof(unsigned int));
		cf11Copy.len = cf11.len;
		for (int i = 0; i < cf11.len; i++)
		{
			cf11Copy.dat[i] = cf11.dat[i];
		}


		//CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);
		//CR_TEST_RET crTest1 = cr_test(&crTest.cf1, thdt, &crTest.nc1);
		//CR_TEST_RET crTest2 = cr_test(&crTest.cf2, thdt, &crTest.nc2);
		//CR_TEST_RET crTest3 = cr_test(&cf12, thdq, &nc12);
		//CR_TEST_RET crTest4 = cr_test(&crTest3.cf1, thdt, &crTest3.nc1);
		//CR_TEST_RET crTest5 = cr_test(&crTest3.cf2, thdt, &crTest3.nc2);
		CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);

		Uint32_Dat crTestCf1Copy;//调用之后数据会被改变，所以需要拷贝
		crTestCf1Copy.len = crTest.cf1.len;
		crTestCf1Copy.dat = calloc(crTest.cf1.len, sizeof(unsigned int));
		for (int i = 0; i < crTest.cf1.len; i++)
		{
			crTestCf1Copy.dat[i] = crTest.cf1.dat[i];
		}
		CR_TEST_RET crTest1 = cr_test(&crTestCf1Copy, thdt, &crTest.nc1);
		//free(crTestCf1Copy.dat);

		Uint32_Dat crTestCf2Copy;//调用之后数据会被改变，所以需要拷贝
		crTestCf2Copy.len = crTest.cf2.len;
		crTestCf2Copy.dat = calloc(crTest.cf2.len, sizeof(unsigned int));
		for (int i = 0; i < crTest.cf2.len; i++)
		{
			crTestCf2Copy.dat[i] = crTest.cf2.dat[i];
		}
		CR_TEST_RET crTest2 = cr_test(&crTestCf2Copy, thdt, &crTest.nc2);
		//free(crTestCf2Copy.dat);

		CR_TEST_RET crTest3 = cr_test(&cf12, thdq, &nc12);

		Uint32_Dat crTest3Cf1Copy;//调用之后数据会被改变，所以需要拷贝
		crTest3Cf1Copy.len = crTest3.cf1.len;
		crTest3Cf1Copy.dat = calloc(crTest3.cf1.len, sizeof(unsigned int));
		for (int i = 0; i < crTest3.cf1.len; i++)
		{
			crTest3Cf1Copy.dat[i] = crTest3.cf1.dat[i];
		}
		CR_TEST_RET crTest4 = cr_test(&crTest3Cf1Copy, thdt, &crTest3.nc1);
		//free(crTest3Cf1Copy.dat);

		Uint32_Dat crTest3Cf2Copy;//调用之后数据会被改变，所以需要拷贝
		crTest3Cf2Copy.len = crTest3.cf2.len;
		crTest3Cf2Copy.dat = calloc(crTest3.cf2.len, sizeof(unsigned int));
		for (int i = 0; i < crTest3.cf2.len; i++)
		{
			crTest3Cf2Copy.dat[i] = crTest3.cf2.dat[i];
		}
		CR_TEST_RET crTest5 = cr_test(&crTest3Cf2Copy, thdt, &crTest3.nc2);
		//free(crTest3Cf2Copy.dat);

		if (crTest.sumcr > 12)
		{
			L2 = 1;

			if (crTest1.sumcr > 12)
			{
				L3 = 1;
				if (crTest2.sumcr > 12)
				{
					L3 = 2;
				}

			}
			if (crTest3.sumcr > 12)
			{
				L2 = 2;

				if ((crTest4.sumcr > 16) && (L3 == 2))
				{
					L3 = 3;

					if (crTest5.sumcr > 16)
					{
						L3 = 4;
					}

				}
			}
		}

		SFC sfc1 = SFcode(L2 + 1, 3);
		binBytePos = ptr & 7;
		rem.a = sfc1.code;
		rem.a = rem.a << (16 - binBytePos - sfc1.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc1.lb;

		SFC sfc2 = SFcode(L3 + 1, 5);
		binBytePos = ptr & 7;
		rem.a = sfc2.code;
		rem.a = rem.a << (16 - binBytePos - sfc2.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc2.lb;

		en_r_cr(&r1, lB1, &nr1, scr1, (float)scrc1, pr, &z1Copy);//binz1

		if (crTest.sumcr > 12)
		{
			//L2 = 1;
			en_r_cr(&(crTest.r), crTest.lB, &(crTest.nr), crTest.sumcr, crTest.sumcrc, crTest.pr, &(crTest.z));//bin_r2A
			if (crTest1.sumcr > 12)
			{
				//L3 = 1;
				en_r_cr(&(crTest1.r), crTest1.lB, &(crTest1.nr), crTest1.sumcr, crTest1.sumcrc, crTest1.pr, &(crTest1.z));//bin_r3A
				en_2n_code_new(&crTest1.cf1, &crTest1.nc1, typ, thdt);//bin31
				en_2n_code_new(&crTest1.cf2, &crTest1.nc2, typ, thdt);//bin31
				if (crTest2.sumcr > 12)
				{
					//L3 = 2;
					en_r_cr(&(crTest2.r), crTest2.lB, &(crTest2.nr), crTest2.sumcr, crTest2.sumcrc, crTest2.pr, &(crTest2.z));//bin_r3B
					en_2n_code_new(&crTest2.cf1, &crTest2.nc1, typ, thdt);//bin33
					en_2n_code_new(&crTest2.cf1, &crTest2.nc2, typ, thdt);//bin33
				}
				else
				{
					en_2n_code_new(&crTest.cf2, &crTest.nc2, typ, thdq);//bin22
				}
			}
			else
			{
				en_2n_code_new(&crTest.cf1, &crTest.nc1, typ, thdq);//bin21
				en_2n_code_new(&crTest.cf2, &crTest.nc2, typ, thdq);//bin22
			}

			if (crTest3.sumcr > 12)
			{
				//L2 = 2;
				en_r_cr(&(crTest3.r), crTest3.lB, &(crTest3.nr), crTest3.sumcr, crTest3.sumcrc, crTest3.pr, &(crTest3.z));//bin_r2B

				if ((crTest4.sumcr > 16) && (L3 == 2))
				{
					//L3 = 3;
					en_r_cr(&(crTest4.r), crTest4.lB, &(crTest4.nr), crTest4.sumcr, crTest4.sumcrc, crTest4.pr, &(crTest4.z));//bin_r3C
					en_2n_code_new(&crTest4.cf1, &crTest4.nc1, typ, thdt);//bin35
					en_2n_code_new(&crTest4.cf2, &crTest4.nc2, typ, thdt);//bin35
					if (crTest5.sumcr > 16)
					{
						//L3 = 4;
						en_r_cr(&(crTest5.r), crTest5.lB, &(crTest5.nr), crTest5.sumcr, crTest5.sumcrc, crTest5.pr, &(crTest5.z));//bin24
						en_2n_code_new(&crTest5.cf1, &crTest5.nc1, typ, thdt);//bin37
						en_2n_code_new(&crTest5.cf2, &crTest5.nc2, typ, thdt);//bin37
					}
					else
					{
						en_2n_code_new(&crTest3.cf2, &crTest3.nc2, typ, thdq);//bin24
					}

				}
				else
				{
					en_2n_code_new(&crTest3.cf1, &crTest3.nc1, typ, thdq);//bin21
					en_2n_code_new(&crTest3.cf2, &crTest3.nc2, typ, thdq);//bin22
				}

			}
			else
			{
				en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin12
			}


		}
		else
		{
			en_2n_code_new(&cf11Copy, &nc11, typ, thdh);//bin11
			en_2n_code_new(&cf12Copy, &nc12, typ, thdh);//bin12
		}


		//CR_TEST_RET crTest = cr_test(&cf11, thdq, &nc11);
		//CR_TEST_RET crTest1 = cr_test(&crTest.cf1, thdt, &crTest.nc1);
		//CR_TEST_RET crTest2 = cr_test(&crTest.cf2, thdt, &crTest.nc2);
		//CR_TEST_RET crTest3 = cr_test(&cf12, thdq, &nc12);
		//CR_TEST_RET crTest4 = cr_test(&crTest3.cf1, thdt, &crTest3.nc1);
		//CR_TEST_RET crTest5 = cr_test(&crTest3.cf2, thdt, &crTest3.nc2);


		free(sepCfk.sep);
		free(sepCfk.rw);
		free(crTest.cf1.dat);
		free(crTest.cf2.dat);
		free(crTest.nc1.dat);
		free(crTest.nc2.dat);
		free(crTest.nr.dat);
		free(crTest.r.dat);
		free(cf11Copy.dat);
		free(cf12Copy.dat);
		free(nc11.dat);
		free(nc12.dat);
		free(z1Copy.dat);


		free(crTest1.cf1.dat);
		free(crTest1.cf2.dat);
		free(crTest1.nc1.dat);
		free(crTest1.nc2.dat);
		free(crTest1.nr.dat);
		free(crTest1.r.dat);

		free(crTest2.cf1.dat);
		free(crTest2.cf2.dat);
		free(crTest2.nc1.dat);
		free(crTest2.nc2.dat);
		free(crTest2.nr.dat);
		free(crTest2.r.dat);

		//free(crTest3.cf1.dat);
		free(crTest3.cf2.dat);
		free(crTest3.nc1.dat);
		free(crTest3.nc2.dat);
		free(crTest3.nr.dat);
		free(crTest3.r.dat);

		free(crTest4.cf1.dat);
		free(crTest4.cf2.dat);
		free(crTest4.nc1.dat);
		free(crTest4.nc2.dat);
		free(crTest4.nr.dat);
		free(crTest4.r.dat);

		free(crTest5.cf1.dat);
		free(crTest5.cf2.dat);
		free(crTest5.nc1.dat);
		free(crTest5.nc2.dat);
		free(crTest5.nr.dat);
		free(crTest5.r.dat);
	}
}