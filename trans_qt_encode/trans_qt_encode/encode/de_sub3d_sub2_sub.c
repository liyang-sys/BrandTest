#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

DE_S_S_SUB de_sub3d_sub2_sub(uchar *biny, int thd, int len, int lenbinbit)
{
	DE_S_S_SUB ret;
	ret.z.dat = NULL;
	ret.z.len = 0;
	ret.cfk.dat = NULL;
	ret.cfk.len = 0;
	extern int ptr;
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
	int thdh = thd / 2, thdq = thd / 4, thdt = thd / 8;
	DEC zDEC = de_r_cr(biny, len, lenbinbit);
	ret.z.dat = zDEC.z;
	ret.z.len = zDEC.lenzbit;

	Uint8_Dat binZ;
	binZ.dat = zDEC.z;
	binZ.len = zDEC.lenzbit;
	int len1 = zDEC.lenzbit - sumBin(&binZ);

	if (thd == 1)
	{
		ret.cfk.len = len1;
		ret.cfk.dat = (unsigned int *)calloc(len1, sizeof(unsigned int));
		for (int i = 0; i < len1; i++)
		{
			ret.cfk.dat[i] = 1;
		}
	}
	else if (thd == 2)
	{
		DEC z1 = de_r_cr(biny, len1, lenbinbit);
		ret.cfk.dat = calloc(z1.lenzbit, sizeof(unsigned int));
		ret.cfk.len = z1.lenzbit;
		int indexZ1 = 0;
		for (int i = 0; i < z1.lenzbit; i++)
		{
			indexZ1 = i & 7;
			if (((z1.z[i >> 3] >> (7 - indexZ1)) & 1))
			{
				ret.cfk.dat[i] = 2;
			}
			else
			{
				ret.cfk.dat[i] = 1;
			}
		}
		free(z1.z);
	}
	else if (thd == 3)
	{
		int thdh = 2;
		binBytePos = ptr & 7;
		int L2 = (biny[ptr >> 3] >> (7 - binBytePos)) & 1;
		ptr++;
		DEC z1 = de_r_cr(biny, len1, lenbinbit);
		Uint8_Dat z;
		z.dat = z1.z;
		z.len = z1.lenzbit;
		int lenw = sumBin(&z);
		int lenk = z1.lenzbit - lenw;
		Uint32_Dat rw;
		rw.len = lenw;
		rw.dat = calloc(rw.len, sizeof(unsigned int));
		for (int i = 0; i < rw.len; i++)
		{
			rw.dat[i] = 1;
		}
		if (L2 == 0)
		{
			typ = 1;
			DEC cf11 = de_2n_code_new(biny, lenk, thdh, typ);
			ret.cfk.dat = separate_inv(z1.z, rw.dat, cf11.r, len1, thdh);
			ret.cfk.len = len1;
			free(cf11.r);
		}
		else
		{
			DEC cf11 = de_r_cr(biny, lenk, lenbinbit);
			cf11.r = (unsigned int*)calloc(cf11.lenzbit, sizeof(unsigned int));
			int indexZ1 = 0;
			for (int i = 0; i < cf11.lenzbit; i++)
			{
				indexZ1 = i & 7;
				if (((cf11.z[i >> 3] >> (7 - indexZ1)) & 1))
				{
					cf11.r[i] = 2;
				}
				else
				{
					cf11.r[i] = 1;
				}
			}

			ret.cfk.dat = separate_inv(z1.z, rw.dat, cf11.r, len1, thdh);
			ret.cfk.len = len1;
			free(cf11.z);
			free(cf11.r);
			//free(cf11.z);
		}
		free(rw.dat);
		free(z1.z);
	}
	else if (thd == 4)
	{
		DES L2 = deSFcode(biny, 3);
		L2.sym--;
		DEC z1 = de_r_cr(biny, len1, lenbinbit);
		Uint8_Dat z;
		z.dat = z1.z;
		z.len = z1.lenzbit;
		int len12 = sumBin(&z);
		int len11 = z1.lenzbit - len12;
		if (L2.sym == 0)
		{
			DEC cf11 = de_2n_code_new(biny, len11, thdh, typ);
			DEC cf12 = de_2n_code_new(biny, len12, thdh, typ);
			ret.cfk.dat = separate_inv(z1.z, cf12.r, cf11.r, len1, thdh);
			ret.cfk.len = len1;
			free(cf11.r);
			free(cf12.r);
		}
		else
		{
			//Œ¥≤‚ ‘
			if (L2.sym == 1)
			{
				DEC zCf11 = de_r_cr(biny, len11, lenbinbit);
				Uint32_Dat cf11;
				cf11.dat = (unsigned int*)calloc(zCf11.lenzbit, sizeof(unsigned int));
				cf11.len = zCf11.lenzbit;
				int indexZ1 = 0;
				for (int i = 0; i < zCf11.lenzbit; i++)
				{
					indexZ1 = i & 7;
					if (((zCf11.z[i >> 3] >> (7 - indexZ1)) & 1))
					{
						cf11.dat[i] = 2;
					}
					else
					{
						cf11.dat[i] = 1;
					}
				}
				DEC cf12 = de_2n_code_new(biny, len12, thdh, typ);
				ret.cfk.dat = separate_inv(z1.z, cf12.r, cf11.dat, len1, thdh);
				ret.cfk.len = len1;
				free(cf11.dat);
				free(cf12.r);
				free(zCf11.z);
			}
			else
			{
				DEC zCf11 = de_r_cr(biny, len11, lenbinbit);
				Uint32_Dat cf11;
				cf11.dat = (unsigned int*)calloc(zCf11.lenzbit, sizeof(unsigned int));
				cf11.len = zCf11.lenzbit;
				int indexZ1 = 0;
				for (int i = 0; i < zCf11.lenzbit; i++)
				{
					indexZ1 = i & 7;
					if (((zCf11.z[i >> 3] >> (7 - indexZ1)) & 1))
					{
						cf11.dat[i] = 2;
					}
					else
					{
						cf11.dat[i] = 1;
					}
				}

				DEC zCf12 = de_r_cr(biny, len12, lenbinbit);
				Uint32_Dat cf12;
				cf12.dat = (unsigned int*)calloc(zCf12.lenzbit, sizeof(unsigned int));
				cf12.len = zCf12.lenzbit;
				indexZ1 = 0;
				for (int i = 0; i < zCf12.lenzbit; i++)
				{
					indexZ1 = i & 7;
					if (((zCf12.z[i >> 3] >> (7 - indexZ1)) & 1))
					{
						cf12.dat[i] = 2;
					}
					else
					{
						cf12.dat[i] = 1;
					}
				}
				ret.cfk.dat = separate_inv(z1.z, cf12.dat, cf11.dat, len1, thdh);
				ret.cfk.len = len1;
				free(cf11.dat);
				free(zCf11.z);
				free(cf12.dat);
				free(zCf12.z);
			}

		}

	}
	else if (thd == 6)
	{
		thdh = 3;
		thdq = 2;
		int typ2 = 1;
		DES L2 = deSFcode(biny, 3);
		L2.sym--;
		DES L3 = deSFcode(biny, 3);
		L3.sym--;
		DEC z1 = de_r_cr(biny, len1, lenbinbit);
		Uint8_Dat z;
		z.dat = z1.z;
		z.len = z1.lenzbit;
		int len12 = sumBin(&z);
		int len11 = z1.lenzbit - len12;

		if (L2.sym == 0)
		{
			DEC cf11 = de_2n_code_new(biny, len11, thdh, typ);
			DEC cf12 = de_2n_code_new(biny, len12, thdh, typ);
			ret.cfk.dat = separate_inv(z1.z, cf12.r, cf11.r, z1.lenzbit, thdh);
			ret.cfk.len = z1.lenzbit;
			free(cf11.r);
			free(cf12.r);
		}
		else
		{
			//Œ¥≤‚ ‘
			if (L2.sym == 1)
			{
				DEC z21 = de_r_cr(biny, len11, lenbinbit);
				Uint8_Dat z21U8v;
				z21U8v.dat = z21.z;
				z21U8v.len = z21.lenzbit;
				int len31 = sumBin(&z21U8v);
				int len3Nd2 = z21.lenzbit - len31;
				if (L3.sym < 1)
				{
					DEC cf21 = de_2n_code_new(biny, len3Nd2, thdq, typ2);
					Uint32_Dat cf11;
					cf11.dat = NULL;
					cf11.len = z21.lenzbit;
					Uint32_Dat temp1;
					temp1.dat = calloc(len31, sizeof(unsigned int));
					temp1.len = len31;
					for (int i = 0; i < len31; i++)
					{
						temp1.dat[i] = 1;
					}
					cf11.dat = separate_inv(z21.z, temp1.dat, cf21.r, z21.lenzbit, thdq);
					DEC cf12 = de_2n_code_new(biny, len12, thdh, typ);

					ret.cfk.dat = separate_inv(z1.z, cf12.r, cf11.dat, z1.lenzbit, thdh);
					ret.cfk.len = z1.lenzbit;

					free(temp1.dat);
					free(cf12.r);
					free(cf11.dat);
					free(cf21.r);
				}
				else
				{
					DEC zCf21 = de_r_cr(biny, len3Nd2, lenbinbit);
					Uint32_Dat cf21;
					cf21.dat = (unsigned int*)calloc(zCf21.lenzbit, sizeof(unsigned int));
					cf21.len = zCf21.lenzbit;
					int indexZ1 = 0;
					for (int i = 0; i < zCf21.lenzbit; i++)
					{
						indexZ1 = i & 7;
						if (((zCf21.z[i >> 3] >> (7 - indexZ1)) & 1))
						{
							cf21.dat[i] = 2;
						}
						else
						{
							cf21.dat[i] = 1;
						}
					}
					Uint32_Dat cf11;
					cf11.dat = NULL;
					cf11.len = z21.lenzbit;
					Uint32_Dat temp1;
					temp1.dat = calloc(len31, sizeof(unsigned int));
					temp1.len = len31;
					for (int i = 0; i < len31; i++)
					{
						temp1.dat[i] = 1;
					}
					cf11.dat = separate_inv(z21.z, temp1.dat, cf21.dat, z21.lenzbit, thdq);
					DEC cf12 = de_2n_code_new(biny, len12, thdh, typ);

					ret.cfk.dat = separate_inv(z1.z, cf12.r, cf11.dat, z1.lenzbit, thdh);
					ret.cfk.len = z1.lenzbit;

					free(temp1.dat);
					free(cf12.r);
					free(cf11.dat);
					free(zCf21.z);
					free(cf21.dat);
				}
				free(z21.z);
			}
			else
			{

				Uint32_Dat cf11;
				cf11.dat = NULL;
				DEC z21 = de_r_cr(biny, len11, lenbinbit);
				Uint8_Dat z21U8v;
				z21U8v.dat = z21.z;
				z21U8v.len = z21.lenzbit;
				int len31 = sumBin(&z21U8v);
				int len3Nd2 = z21.lenzbit - len31;
				if (L3.sym < 1)
				{
					DEC cf21 = de_2n_code_new(biny, len3Nd2, thdq, typ2);
					Uint32_Dat temp1;
					temp1.dat = calloc(len31, sizeof(unsigned int));
					temp1.len = len31;
					for (int i = 0; i < len31; i++)
					{
						temp1.dat[i] = 1;
					}

					cf11.len = z21.lenzbit;
					cf11.dat = separate_inv(z21.z, temp1.dat, cf21.r, z21.lenzbit, thdq);
					free(temp1.dat);
					free(cf21.r);
				}
				else
				{
					DEC zCf21 = de_r_cr(biny, len3Nd2, lenbinbit);
					Uint32_Dat cf21;
					cf21.dat = (unsigned int*)calloc(zCf21.lenzbit, sizeof(unsigned int));
					cf21.len = zCf21.lenzbit;
					int indexZ1 = 0;
					for (int i = 0; i < zCf21.lenzbit; i++)
					{
						indexZ1 = i & 7;
						if (((zCf21.z[i >> 3] >> (7 - indexZ1)) & 1))
						{
							cf21.dat[i] = 2;
						}
						else
						{
							cf21.dat[i] = 1;
						}
					}

					Uint32_Dat temp1;
					temp1.dat = calloc(len31, sizeof(unsigned int));
					temp1.len = len31;
					for (int i = 0; i < len31; i++)
					{
						temp1.dat[i] = 1;
					}
					cf11.len = z21.lenzbit;
					cf11.dat = separate_inv(z21.z, temp1.dat, cf21.dat, z21.lenzbit, thdq);
					free(temp1.dat);
					free(zCf21.z);
					free(cf21.dat);
				}

				Uint32_Dat cf12;
				cf12.dat = NULL;
				DEC z22 = de_r_cr(biny, len12, lenbinbit);
				Uint8_Dat z22U8v;
				z22U8v.dat = z22.z;
				z22U8v.len = z22.lenzbit;
				int len32 = sumBin(&z22U8v);
				int len32Nd2 = z22.lenzbit - len32;
				if (L3.sym < 2)
				{
					DEC cf22 = de_2n_code_new(biny, len32Nd2, thdq, typ2);
					cf12.len = z22.lenzbit;
					Uint32_Dat temp1;
					temp1.dat = calloc(len31, sizeof(unsigned int));
					temp1.len = len31;
					for (int i = 0; i < len31; i++)
					{
						temp1.dat[i] = 1;
					}
					//
					cf12.dat = separate_inv(z22.z, temp1.dat, cf22.r, z22.lenzbit, thdq);
					free(temp1.dat);
					free(cf22.r);
				}
				else
				{
					DEC zCf22 = de_r_cr(biny, len32Nd2, lenbinbit);
					Uint32_Dat cf22;
					cf22.dat = (unsigned int*)calloc(zCf22.lenzbit, sizeof(unsigned int));
					cf22.len = zCf22.lenzbit;
					int indexZ1 = 0;
					for (int i = 0; i < zCf22.lenzbit; i++)
					{
						indexZ1 = i & 7;
						if (((zCf22.z[i >> 3] >> (7 - indexZ1)) & 1))
						{
							cf22.dat[i] = 2;
						}
						else
						{
							cf22.dat[i] = 1;
						}
					}
					Uint32_Dat temp1;
					temp1.dat = calloc(len31, sizeof(unsigned int));
					temp1.len = len31;
					for (int i = 0; i < len31; i++)
					{
						temp1.dat[i] = 1;
					}
					cf12.len = z22.lenzbit;
					cf12.dat = separate_inv(z22.z, temp1.dat, cf22.dat, z22.lenzbit, thdq);
					free(temp1.dat);
					free(zCf22.z);
					free(cf22.dat);
				}
				ret.cfk.dat = separate_inv(z1.z, cf12.dat, cf11.dat, z1.lenzbit, thdh);
				ret.cfk.len = z1.lenzbit;
				free(z21.z);
				free(z22.z);
				free(cf11.dat);
				free(cf12.dat);
			}
		}
		free(z1.z);
	}
	else if (thd == 8)
	{
		DES L2 = deSFcode(biny, 3);
		L2.sym--;
		DES L3 = deSFcode(biny, 5);
		L3.sym--;
		DEC z1 = de_r_cr(biny, len1, lenbinbit);
		Uint8_Dat z;
		z.dat = z1.z;
		z.len = z1.lenzbit;
		int len12 = sumBin(&z);
		int len11 = z1.lenzbit - len12;

		if (L2.sym == 0)
		{
			DEC cf11 = de_2n_code_new(biny, len11, thdh, typ);
			DEC cf12 = de_2n_code_new(biny, len12, thdh, typ);
			ret.cfk.dat = separate_inv(z1.z, cf12.r, cf11.r, z1.lenzbit, thdh);
			ret.cfk.len = z1.lenzbit;
			free(cf11.r);
			free(cf12.r);
		}
		else
		{
			if (L2.sym == 1)
			{
				int na = 1, nb = 2;
				DEC z21 = de_r_cr(biny, len11, lenbinbit);
				Uint8_Dat z21U8v;
				z21U8v.dat = z21.z;
				z21U8v.len = z21.lenzbit;
				int len32 = sumBin(&z21U8v);
				int len31 = z21.lenzbit - len32;
				Uint32_Dat cf21;
				cf21.dat = NULL;
				cf21.len = 0;

				Uint32_Dat cf11;
				cf11.dat = NULL;
				cf11.len = 0;
				//nd3 == 1
				if (L3.sym < 1)
				{
					DEC cf21Dec = de_2n_code_new(biny, len31, thdq, typ);
					cf21.dat = cf21Dec.r;
					cf21.len = cf21Dec.lenr;
				}
				else
				{
					DEC zCf31 = de_r_cr(biny, len31, lenbinbit);
					cf21.dat = (unsigned int*)calloc(zCf31.lenzbit, sizeof(unsigned int));
					cf21.len = zCf31.lenzbit;
					int indexZ1 = 0;
					for (int i = 0; i < zCf31.lenzbit; i++)
					{
						indexZ1 = i & 7;
						if (((zCf31.z[i >> 3] >> (7 - indexZ1)) & 1))
						{
							cf21.dat[i] = 2;
						}
						else
						{
							cf21.dat[i] = 1;
						}
					}
					free(zCf31.z);
				}
				//nd3 == 2
				Uint32_Dat cf22;
				cf22.dat = NULL;
				cf22.len = 0;
				if (L3.sym < 2)
				{
					DEC cf22Dec = de_2n_code_new(biny, len32, thdq, typ);
					cf22.dat = cf22Dec.r;
					cf22.len = cf22Dec.lenr;
				}
				else
				{
					DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					cf22.dat = (unsigned int*)calloc(zCf32.lenzbit, sizeof(unsigned int));
					cf22.len = zCf32.lenzbit;
					int indexZ1 = 0;
					for (int i = 0; i < zCf32.lenzbit; i++)
					{
						indexZ1 = i & 7;
						if (((zCf32.z[i >> 3] >> (7 - indexZ1)) & 1))
						{
							cf22.dat[i] = 2;
						}
						else
						{
							cf22.dat[i] = 1;
						}
					}
					free(zCf32.z);
				}

				cf11.len = z21.lenzbit;
				cf11.dat = separate_inv(z21.z, cf22.dat, cf21.dat, z21.lenzbit, thdq);

				Uint32_Dat cf12;
				cf12.dat = NULL;
				cf12.len = 0;

				DEC cf12Dec = de_2n_code_new(biny, len12, thdh, typ);
				cf12.dat = cf12Dec.r;
				cf12.len = cf12Dec.lenr;

				ret.cfk.dat = separate_inv(z1.z, cf12.dat, cf11.dat, z1.lenzbit, thdh);
				ret.cfk.len = z1.lenzbit;

				free(z21.z);
				if (cf21.dat != NULL)
				{
					free(cf21.dat);
				}
				if (cf22.dat != NULL)
				{
					free(cf22.dat);
				}
				if (cf11.dat != NULL)
				{
					free(cf11.dat);
				}
				if (cf12.dat != NULL)
				{
					free(cf12.dat);
				}
			}
			else
			{
				int na = 1, nb = 2;
				DEC z21 = de_r_cr(biny, len11, lenbinbit);
				Uint8_Dat z21U8v;
				z21U8v.dat = z21.z;
				z21U8v.len = z21.lenzbit;
				int len32 = sumBin(&z21U8v);
				int len31 = z21.lenzbit - len32;
				Uint32_Dat cf21;
				cf21.dat = NULL;
				cf21.len = 0;

				Uint32_Dat cf11;
				cf11.dat = NULL;
				cf11.len = 0;
				//nd3 == 1
				if (L3.sym < 1)
				{
					DEC cf21Dec = de_2n_code_new(biny, len31, thdq, typ);
					cf21.dat = cf21Dec.r;
					cf21.len = cf21Dec.lenr;
				}
				else
				{
					DEC zCf31 = de_r_cr(biny, len31, lenbinbit);
					cf21.dat = (unsigned int*)calloc(zCf31.lenzbit, sizeof(unsigned int));
					cf21.len = zCf31.lenzbit;
					int indexZ1 = 0;
					for (int i = 0; i < zCf31.lenzbit; i++)
					{
						indexZ1 = i & 7;
						if (((zCf31.z[i >> 3] >> (7 - indexZ1)) & 1))
						{
							cf21.dat[i] = 2;
						}
						else
						{
							cf21.dat[i] = 1;
						}
					}
					free(zCf31.z);
				}
				//nd3 == 2
				Uint32_Dat cf22;
				cf22.dat = NULL;
				cf22.len = 0;
				if (L3.sym < 2)
				{
					DEC cf22Dec = de_2n_code_new(biny, len32, thdq, typ);
					cf22.dat = cf22Dec.r;
					cf22.len = cf22Dec.lenr;
				}
				else
				{
					DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					cf22.dat = (unsigned int*)calloc(zCf32.lenzbit, sizeof(unsigned int));
					cf22.len = zCf32.lenzbit;
					int indexZ1 = 0;
					for (int i = 0; i < zCf32.lenzbit; i++)
					{
						indexZ1 = i & 7;
						if (((zCf32.z[i >> 3] >> (7 - indexZ1)) & 1))
						{
							cf22.dat[i] = 2;
						}
						else
						{
							cf22.dat[i] = 1;
						}
					}
					free(zCf32.z);
				}

				cf11.len = z21.lenzbit;
				cf11.dat = separate_inv(z21.z, cf22.dat, cf21.dat, z21.lenzbit, thdq);



				na = 3, nb = 4;
				DEC z22 = de_r_cr(biny, len12, lenbinbit);
				Uint8_Dat z22U8v;
				z22U8v.dat = z22.z;
				z22U8v.len = z22.lenzbit;
				int len34 = sumBin(&z22U8v);
				int len33 = z22.lenzbit - len34;

				//nd3 == 3
				Uint32_Dat cf23;
				cf23.dat = NULL;
				cf23.len = 0;

				Uint32_Dat cf12;
				cf12.dat = NULL;
				cf12.len = 0;

				if (L3.sym < 3)
				{
					DEC cf22Dec = de_2n_code_new(biny, len33, thdq, typ);
					cf23.dat = cf22Dec.r;
					cf23.len = cf22Dec.lenr;
				}
				else
				{
					DEC zCf33 = de_r_cr(biny, len33, lenbinbit);
					cf23.dat = (unsigned int*)calloc(zCf33.lenzbit, sizeof(unsigned int));
					cf23.len = zCf33.lenzbit;
					int indexZ1 = 0;
					for (int i = 0; i < zCf33.lenzbit; i++)
					{
						indexZ1 = i & 7;
						if (((zCf33.z[i >> 3] >> (7 - indexZ1)) & 1))
						{
							cf23.dat[i] = 2;
						}
						else
						{
							cf23.dat[i] = 1;
						}
					}
					free(zCf33.z);
				}

				//nd3 == 4
				Uint32_Dat cf24;
				cf24.dat = NULL;
				cf24.len = 0;

				if (L3.sym < 4)
				{
					DEC cf24Dec = de_2n_code_new(biny, len34, thdq, typ);
					cf24.dat = cf24Dec.r;
					cf24.len = cf24Dec.lenr;
				}
				else
				{
					DEC zCf34 = de_r_cr(biny, len34, lenbinbit);
					cf24.dat = (unsigned int*)calloc(zCf34.lenzbit, sizeof(unsigned int));
					cf24.len = zCf34.lenzbit;
					int indexZ1 = 0;
					for (int i = 0; i < zCf34.lenzbit; i++)
					{
						indexZ1 = i & 7;
						if (((zCf34.z[i >> 3] >> (7 - indexZ1)) & 1))
						{
							cf24.dat[i] = 2;
						}
						else
						{
							cf24.dat[i] = 1;
						}
					}
					free(zCf34.z);
				}

				cf12.len = z22.lenzbit;
				cf12.dat = separate_inv(z22.z, cf24.dat, cf23.dat, z22.lenzbit, thdq);

				ret.cfk.dat = separate_inv(z1.z, cf12.dat, cf11.dat, z1.lenzbit, thdh);
				ret.cfk.len = z1.lenzbit;

				free(z21.z);
				if (cf21.dat != NULL)
				{
					free(cf21.dat);
				}
				if (cf22.dat != NULL)
				{
					free(cf22.dat);
				}
				if (cf11.dat != NULL)
				{
					free(cf11.dat);
				}
				free(z22.z);
				free(cf23.dat);
				free(cf24.dat);
				free(cf12.dat);
			}
		}
	}
	else if (thd == 12)
	{
		thdh = 6;
		thdq = 3;
		int thdt = 2;
		int typ3 = 1;

		DES L2 = deSFcode(biny, 3);
		L2.sym--;
		DES L3 = deSFcode(biny, 5);
		L3.sym--;
		DEC z1 = de_r_cr(biny, len1, lenbinbit);
		Uint8_Dat z;
		z.dat = z1.z;
		z.len = z1.lenzbit;
		int len12 = sumBin(&z);
		int len11 = z1.lenzbit - len12;

		if (L2.sym == 0)
		{
			DEC cf11 = de_2n_code_new(biny, len11, thdh, typ);
			DEC cf12 = de_2n_code_new(biny, len12, thdh, typ);
			ret.cfk.dat = separate_inv(z1.z, cf12.r, cf11.r, z1.lenzbit, thdh);
			ret.cfk.len = z1.lenzbit;
			free(cf11.r);
			free(cf12.r);
		}
		else
		{
			if (L2.sym == 1)
			{
				int na = 1, nb = 2;
				DEC z21 = de_r_cr(biny, len11, lenbinbit);
				Uint8_Dat z21U8v;
				z21U8v.dat = z21.z;
				z21U8v.len = z21.lenzbit;
				int len32 = sumBin(&z21U8v);
				int len31 = z21.lenzbit - len32;
				Uint32_Dat cf21;
				cf21.dat = NULL;
				cf21.len = 0;

				Uint32_Dat cf11;
				cf11.dat = NULL;
				cf11.len = 0;
				//nd3 == 1
				if (L3.sym < 1)
				{
					DEC cf21Dec = de_2n_code_new(biny, len31, thdq, typ);
					cf21.dat = cf21Dec.r;
					cf21.len = cf21Dec.lenr;
				}
				else
				{
					DEC zCf31 = de_r_cr(biny, len31, lenbinbit);
					Uint8_Dat z31U8v;
					z31U8v.dat = zCf31.z;
					z31U8v.len = zCf31.lenzbit;
					int lenw3 = sumBin(&z31U8v);
					int lenk3 = zCf31.lenzbit - lenw3;
					DEC cf3k1Dec = de_2n_code_new(biny, lenk3, thdt, typ3);
					Uint32_Dat cf3w1;
					cf3w1.len = lenw3;
					cf3w1.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					for (int i = 0; i < lenw3; i++)
					{
						cf3w1.dat[i] = 1;
					}
					cf21.len = zCf31.lenzbit;
					cf21.dat = separate_inv(zCf31.z, cf3w1.dat, cf3k1Dec.r, zCf31.lenzbit, thdt);
					free(zCf31.z);
					free(cf3k1Dec.r);
					free(cf3w1.dat);
				}
				//nd3 == 2
				Uint32_Dat cf22;
				cf22.dat = NULL;
				cf22.len = 0;
				if (L3.sym < 2)
				{
					DEC cf22Dec = de_2n_code_new(biny, len32, thdq, typ);
					cf22.dat = cf22Dec.r;
					cf22.len = cf22Dec.lenr;
				}
				else
				{
					//DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					//cf22.dat = (unsigned int*)calloc(zCf32.lenzbit, sizeof(unsigned int));
					//cf22.len = zCf32.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf32.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf32.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf22.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf22.dat[i] = 1;
					//	}
					//}
					//free(zCf32.z);
					DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					Uint8_Dat z32U8v;
					z32U8v.dat = zCf32.z;
					z32U8v.len = zCf32.lenzbit;
					int lenw3 = sumBin(&z32U8v);
					int lenk3 = zCf32.lenzbit - lenw3;
					DEC cf3k2Dec = de_2n_code_new(biny, lenk3, thdt, typ3);
					Uint32_Dat cf3w2;
					cf3w2.len = lenw3;
					cf3w2.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					for (int i = 0; i < lenw3; i++)
					{
						cf3w2.dat[i] = 1;
					}
					cf22.len = zCf32.lenzbit;
					cf22.dat = separate_inv(zCf32.z, cf3w2.dat, cf3k2Dec.r, zCf32.lenzbit, thdt);
					free(zCf32.z);
					free(cf3k2Dec.r);
					free(cf3w2.dat);
				}

				cf11.len = z21.lenzbit;
				cf11.dat = separate_inv(z21.z, cf22.dat, cf21.dat, z21.lenzbit, thdq);

				Uint32_Dat cf12;
				cf12.dat = NULL;
				cf12.len = 0;

				DEC cf12Dec = de_2n_code_new(biny, len12, thdh, typ);
				cf12.dat = cf12Dec.r;
				cf12.len = cf12Dec.lenr;

				ret.cfk.dat = separate_inv(z1.z, cf12.dat, cf11.dat, z1.lenzbit, thdh);
				ret.cfk.len = z1.lenzbit;

				free(z21.z);
				if (cf21.dat != NULL)
				{

					free(cf21.dat);
				}
				if (cf22.dat != NULL)
				{
					free(cf22.dat);
				}
				if (cf11.dat != NULL)
				{
					free(cf11.dat);
				}
				if (cf12.dat != NULL)
				{
					free(cf12.dat);
				}
			}
			else
			{
				int na = 1, nb = 2;
				DEC z21 = de_r_cr(biny, len11, lenbinbit);
				Uint8_Dat z21U8v;
				z21U8v.dat = z21.z;
				z21U8v.len = z21.lenzbit;
				int len32 = sumBin(&z21U8v);
				int len31 = z21.lenzbit - len32;
				Uint32_Dat cf21;
				cf21.dat = NULL;
				cf21.len = 0;

				Uint32_Dat cf11;
				cf11.dat = NULL;
				cf11.len = 0;
				//nd3 == 1
				if (L3.sym < 1)
				{
					DEC cf21Dec = de_2n_code_new(biny, len31, thdq, typ);
					cf21.dat = cf21Dec.r;
					cf21.len = cf21Dec.lenr;
				}
				else
				{
					//DEC zCf31 = de_r_cr(biny, len31, lenbinbit);
					//cf21.dat = (unsigned int*)calloc(zCf31.lenzbit, sizeof(unsigned int));
					//cf21.len = zCf31.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf31.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf31.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf21.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf21.dat[i] = 1;
					//	}
					//}
					//free(zCf31.z);
					DEC zCf31 = de_r_cr(biny, len31, lenbinbit);
					Uint8_Dat z31U8v;
					z31U8v.dat = zCf31.z;
					z31U8v.len = zCf31.lenzbit;
					int lenw3 = sumBin(&z31U8v);
					int lenk3 = zCf31.lenzbit - lenw3;
					DEC cf3k1Dec = de_2n_code_new(biny, lenk3, thdt, typ3);
					Uint32_Dat cf3w1;
					cf3w1.len = lenw3;
					cf3w1.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					for (int i = 0; i < lenw3; i++)
					{
						cf3w1.dat[i] = 1;
					}
					cf21.len = zCf31.lenzbit;
					cf21.dat = separate_inv(zCf31.z, cf3w1.dat, cf3k1Dec.r, zCf31.lenzbit, thdt);
					free(zCf31.z);
					free(cf3k1Dec.r);
					free(cf3w1.dat);
				}
				//nd3 == 2
				Uint32_Dat cf22;
				cf22.dat = NULL;
				cf22.len = 0;
				if (L3.sym < 2)
				{
					DEC cf22Dec = de_2n_code_new(biny, len32, thdq, typ);
					cf22.dat = cf22Dec.r;
					cf22.len = cf22Dec.lenr;
				}
				else
				{
					//DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					//cf22.dat = (unsigned int*)calloc(zCf32.lenzbit, sizeof(unsigned int));
					//cf22.len = zCf32.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf32.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf32.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf22.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf22.dat[i] = 1;
					//	}
					//}
					//free(zCf32.z);
					DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					Uint8_Dat z32U8v;
					z32U8v.dat = zCf32.z;
					z32U8v.len = zCf32.lenzbit;
					int lenw3 = sumBin(&z32U8v);
					int lenk3 = zCf32.lenzbit - lenw3;
					DEC cf3k2Dec = de_2n_code_new(biny, lenk3, thdt, typ3);
					Uint32_Dat cf3w2;
					cf3w2.len = lenw3;
					cf3w2.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					for (int i = 0; i < lenw3; i++)
					{
						cf3w2.dat[i] = 1;
					}
					cf22.len = zCf32.lenzbit;
					cf22.dat = separate_inv(zCf32.z, cf3w2.dat, cf3k2Dec.r, zCf32.lenzbit, thdt);
					free(zCf32.z);
					free(cf3k2Dec.r);
					free(cf3w2.dat);
				}

				cf11.len = z21.lenzbit;
				cf11.dat = separate_inv(z21.z, cf22.dat, cf21.dat, z21.lenzbit, thdq);



				na = 3, nb = 4;
				DEC z22 = de_r_cr(biny, len12, lenbinbit);
				Uint8_Dat z22U8v;
				z22U8v.dat = z22.z;
				z22U8v.len = z22.lenzbit;
				int len34 = sumBin(&z22U8v);
				int len33 = z22.lenzbit - len34;

				//nd3 == 3
				Uint32_Dat cf23;
				cf23.dat = NULL;
				cf23.len = 0;

				Uint32_Dat cf12;
				cf12.dat = NULL;
				cf12.len = 0;

				if (L3.sym < 3)
				{
					DEC cf22Dec = de_2n_code_new(biny, len33, thdq, typ);
					cf23.dat = cf22Dec.r;
					cf23.len = cf22Dec.lenr;
				}
				else
				{
					//DEC zCf33 = de_r_cr(biny, len33, lenbinbit);
					//cf23.dat = (unsigned int*)calloc(zCf33.lenzbit, sizeof(unsigned int));
					//cf23.len = zCf33.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf33.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf33.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf23.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf23.dat[i] = 1;
					//	}
					//}
					//free(zCf33.z);
					DEC zCf33 = de_r_cr(biny, len33, lenbinbit);
					Uint8_Dat z33U8v;
					z33U8v.dat = zCf33.z;
					z33U8v.len = zCf33.lenzbit;
					int lenw3 = sumBin(&z33U8v);
					int lenk3 = zCf33.lenzbit - lenw3;
					DEC cf3k3Dec = de_2n_code_new(biny, lenk3, thdt, typ3);
					Uint32_Dat cf3w3;
					cf3w3.len = lenw3;
					cf3w3.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					for (int i = 0; i < lenw3; i++)
					{
						cf3w3.dat[i] = 1;
					}
					cf23.len = zCf33.lenzbit;
					cf23.dat = separate_inv(zCf33.z, cf3w3.dat, cf3k3Dec.r, zCf33.lenzbit, thdt);
					free(zCf33.z);
					free(cf3k3Dec.r);
					free(cf3w3.dat);
				}

				//nd3 == 4
				Uint32_Dat cf24;
				cf24.dat = NULL;
				cf24.len = 0;
				if (L3.sym < 4)
				{
					DEC cf22Dec = de_2n_code_new(biny, len34, thdq, typ);
					cf24.dat = cf22Dec.r;
					cf24.len = cf22Dec.lenr;
				}
				else
				{
					//DEC zCf33 = de_r_cr(biny, len33, lenbinbit);
					//cf23.dat = (unsigned int*)calloc(zCf33.lenzbit, sizeof(unsigned int));
					//cf23.len = zCf33.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf33.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf33.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf23.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf23.dat[i] = 1;
					//	}
					//}
					//free(zCf33.z);
					DEC zCf34 = de_r_cr(biny, len34, lenbinbit);
					Uint8_Dat z34U8v;
					z34U8v.dat = zCf34.z;
					z34U8v.len = zCf34.lenzbit;
					int lenw3 = sumBin(&z34U8v);
					int lenk3 = zCf34.lenzbit - lenw3;
					DEC cf3k4Dec = de_2n_code_new(biny, lenk3, thdt, typ3);
					Uint32_Dat cf3w4;
					cf3w4.len = lenw3;
					cf3w4.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					for (int i = 0; i < lenw3; i++)
					{
						cf3w4.dat[i] = 1;
					}
					cf24.len = zCf34.lenzbit;
					cf24.dat = separate_inv(zCf34.z, cf3w4.dat, cf3k4Dec.r, zCf34.lenzbit, thdt);
					free(zCf34.z);
					free(cf3k4Dec.r);
					free(cf3w4.dat);
				}




				cf12.len = z22.lenzbit;
				cf12.dat = separate_inv(z22.z, cf24.dat, cf23.dat, z22.lenzbit, thdq);

				ret.cfk.dat = separate_inv(z1.z, cf12.dat, cf11.dat, z1.lenzbit, thdh);
				ret.cfk.len = z1.lenzbit;

				free(z21.z);
				if (cf21.dat != NULL)
				{
					free(cf21.dat);
				}
				if (cf22.dat != NULL)
				{
					free(cf22.dat);
				}
				if (cf11.dat != NULL)
				{
					free(cf11.dat);
				}
				free(z22.z);
				free(cf23.dat);
				free(cf24.dat);
				free(cf12.dat);
			}
		}
	}
	else if (thd >= 16)
	{
		DES L2 = deSFcode(biny, 3);
		L2.sym--;
		DES L3 = deSFcode(biny, 5);
		L3.sym--;
		DEC z1 = de_r_cr(biny, len1, lenbinbit);
		Uint8_Dat z;
		z.dat = z1.z;
		z.len = z1.lenzbit;
		int len12 = sumBin(&z);
		int len11 = z1.lenzbit - len12;

		if (L2.sym == 0)
		{
			DEC cf11 = de_2n_code_new(biny, len11, thdh, typ);
			DEC cf12 = de_2n_code_new(biny, len12, thdh, typ);
			ret.cfk.dat = separate_inv(z1.z, cf12.r, cf11.r, z1.lenzbit, thdh);
			ret.cfk.len = z1.lenzbit;
			free(cf11.r);
			free(cf12.r);
		}
		else
		{
			if (L2.sym == 1)
			{
				int na = 1, nb = 2;
				DEC z21 = de_r_cr(biny, len11, lenbinbit);
				Uint8_Dat z21U8v;
				z21U8v.dat = z21.z;
				z21U8v.len = z21.lenzbit;
				int len32 = sumBin(&z21U8v);
				int len31 = z21.lenzbit - len32;
				Uint32_Dat cf21;
				cf21.dat = NULL;
				cf21.len = 0;

				Uint32_Dat cf11;
				cf11.dat = NULL;
				cf11.len = 0;
				//nd3 == 1
				if (L3.sym < 1)
				{
					DEC cf21Dec = de_2n_code_new(biny, len31, thdq, typ);
					cf21.dat = cf21Dec.r;
					cf21.len = cf21Dec.lenr;
				}
				else
				{
					DEC zCf31 = de_r_cr(biny, len31, lenbinbit);
					Uint8_Dat z31U8v;
					z31U8v.dat = zCf31.z;
					z31U8v.len = zCf31.lenzbit;
					int lenw3 = sumBin(&z31U8v);
					int lenk3 = zCf31.lenzbit - lenw3;
					DEC cf3k1Dec = de_2n_code_new(biny, lenk3, thdt, typ);
					DEC cf3w1Dec = de_2n_code_new(biny, lenw3, thdt, typ);

					cf21.len = zCf31.lenzbit;
					cf21.dat = separate_inv(zCf31.z, cf3w1Dec.r, cf3k1Dec.r, zCf31.lenzbit, thdt);
					free(zCf31.z);
					free(cf3k1Dec.r);
					free(cf3w1Dec.r);
				}
				//nd3 == 2
				Uint32_Dat cf22;
				cf22.dat = NULL;
				cf22.len = 0;
				if (L3.sym < 2)
				{
					DEC cf22Dec = de_2n_code_new(biny, len32, thdq, typ);
					cf22.dat = cf22Dec.r;
					cf22.len = cf22Dec.lenr;
				}
				else
				{
					//DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					//cf22.dat = (unsigned int*)calloc(zCf32.lenzbit, sizeof(unsigned int));
					//cf22.len = zCf32.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf32.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf32.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf22.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf22.dat[i] = 1;
					//	}
					//}
					//free(zCf32.z);
					DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					Uint8_Dat z32U8v;
					z32U8v.dat = zCf32.z;
					z32U8v.len = zCf32.lenzbit;
					int lenw3 = sumBin(&z32U8v);
					int lenk3 = zCf32.lenzbit - lenw3;
					DEC cf3k2Dec = de_2n_code_new(biny, lenk3, thdt, typ);
					//Uint32_Dat cf3w2;
					//cf3w2.len = lenw3;
					//cf3w2.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					//for (int i = 0; i < lenw3; i++)
					//{
					//	cf3w2.dat[i] = 1;
					//}
					DEC cf3w2Dec = de_2n_code_new(biny, lenw3, thdt, typ);
					cf22.len = zCf32.lenzbit;
					cf22.dat = separate_inv(zCf32.z, cf3w2Dec.r, cf3k2Dec.r, zCf32.lenzbit, thdt);
					free(zCf32.z);
					free(cf3k2Dec.r);
					free(cf3w2Dec.r);
				}

				cf11.len = z21.lenzbit;
				cf11.dat = separate_inv(z21.z, cf22.dat, cf21.dat, z21.lenzbit, thdq);

				Uint32_Dat cf12;
				cf12.dat = NULL;
				cf12.len = 0;

				DEC cf12Dec = de_2n_code_new(biny, len12, thdh, typ);
				cf12.dat = cf12Dec.r;
				cf12.len = cf12Dec.lenr;

				ret.cfk.dat = separate_inv(z1.z, cf12.dat, cf11.dat, z1.lenzbit, thdh);
				ret.cfk.len = z1.lenzbit;

				free(z21.z);
				if (cf21.dat != NULL)
				{
					free(cf21.dat);
				}
				if (cf22.dat != NULL)
				{
					free(cf22.dat);
				}
				if (cf11.dat != NULL)
				{
					free(cf11.dat);
				}
				if (cf12.dat != NULL)
				{
					free(cf12.dat);
				}
			}
			else
			{
				int na = 1, nb = 2;
				DEC z21 = de_r_cr(biny, len11, lenbinbit);
				Uint8_Dat z21U8v;
				z21U8v.dat = z21.z;
				z21U8v.len = z21.lenzbit;
				int len32 = sumBin(&z21U8v);
				int len31 = z21.lenzbit - len32;
				Uint32_Dat cf21;
				cf21.dat = NULL;
				cf21.len = 0;

				Uint32_Dat cf11;
				cf11.dat = NULL;
				cf11.len = 0;
				//nd3 == 1
				if (L3.sym < 1)
				{
					DEC cf21Dec = de_2n_code_new(biny, len31, thdq, typ);
					cf21.dat = cf21Dec.r;
					cf21.len = cf21Dec.lenr;
				}
				else
				{
					//DEC zCf31 = de_r_cr(biny, len31, lenbinbit);
					//cf21.dat = (unsigned int*)calloc(zCf31.lenzbit, sizeof(unsigned int));
					//cf21.len = zCf31.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf31.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf31.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf21.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf21.dat[i] = 1;
					//	}
					//}
					//free(zCf31.z);
					DEC zCf31 = de_r_cr(biny, len31, lenbinbit);
					Uint8_Dat z31U8v;
					z31U8v.dat = zCf31.z;
					z31U8v.len = zCf31.lenzbit;
					int lenw3 = sumBin(&z31U8v);
					int lenk3 = zCf31.lenzbit - lenw3;
					DEC cf3k1Dec = de_2n_code_new(biny, lenk3, thdt, typ);
					//Uint32_Dat cf3w1;
					//cf3w1.len = lenw3;
					//cf3w1.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					//for (int i = 0; i < lenw3; i++)
					//{
					//	cf3w1.dat[i] = 1;
					//}
					DEC cf3w1Dec = de_2n_code_new(biny, lenw3, thdt, typ);
					cf21.len = zCf31.lenzbit;
					cf21.dat = separate_inv(zCf31.z, cf3w1Dec.r, cf3k1Dec.r, zCf31.lenzbit, thdt);
					free(zCf31.z);
					free(cf3k1Dec.r);
					free(cf3w1Dec.r);
				}
				//nd3 == 2
				Uint32_Dat cf22;
				cf22.dat = NULL;
				cf22.len = 0;
				if (L3.sym < 2)
				{
					DEC cf22Dec = de_2n_code_new(biny, len32, thdq, typ);
					cf22.dat = cf22Dec.r;
					cf22.len = cf22Dec.lenr;
				}
				else
				{
					//DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					//cf22.dat = (unsigned int*)calloc(zCf32.lenzbit, sizeof(unsigned int));
					//cf22.len = zCf32.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf32.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf32.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf22.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf22.dat[i] = 1;
					//	}
					//}
					//free(zCf32.z);
					DEC zCf32 = de_r_cr(biny, len32, lenbinbit);
					Uint8_Dat z32U8v;
					z32U8v.dat = zCf32.z;
					z32U8v.len = zCf32.lenzbit;
					int lenw3 = sumBin(&z32U8v);
					int lenk3 = zCf32.lenzbit - lenw3;
					DEC cf3k2Dec = de_2n_code_new(biny, lenk3, thdt, typ);
					//Uint32_Dat cf3w2;
					//cf3w2.len = lenw3;
					//cf3w2.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					//for (int i = 0; i < lenw3; i++)
					//{
					//	cf3w2.dat[i] = 1;
					//}
					DEC cf3w2Dec = de_2n_code_new(biny, lenw3, thdt, typ);
					cf22.len = zCf32.lenzbit;
					cf22.dat = separate_inv(zCf32.z, cf3w2Dec.r, cf3k2Dec.r, zCf32.lenzbit, thdt);
					free(zCf32.z);
					free(cf3k2Dec.r);
					free(cf3w2Dec.r);
				}

				cf11.len = z21.lenzbit;
				cf11.dat = separate_inv(z21.z, cf22.dat, cf21.dat, z21.lenzbit, thdq);



				na = 3, nb = 4;
				DEC z22 = de_r_cr(biny, len12, lenbinbit);
				Uint8_Dat z22U8v;
				z22U8v.dat = z22.z;
				z22U8v.len = z22.lenzbit;
				int len34 = sumBin(&z22U8v);
				int len33 = z22.lenzbit - len34;

				//nd3 == 3
				Uint32_Dat cf23;
				cf23.dat = NULL;
				cf23.len = 0;

				Uint32_Dat cf12;
				cf12.dat = NULL;
				cf12.len = 0;

				if (L3.sym < 3)
				{
					DEC cf22Dec = de_2n_code_new(biny, len33, thdq, typ);
					cf23.dat = cf22Dec.r;
					cf23.len = cf22Dec.lenr;
				}
				else
				{
					//DEC zCf33 = de_r_cr(biny, len33, lenbinbit);
					//cf23.dat = (unsigned int*)calloc(zCf33.lenzbit, sizeof(unsigned int));
					//cf23.len = zCf33.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf33.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf33.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf23.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf23.dat[i] = 1;
					//	}
					//}
					//free(zCf33.z);
					DEC zCf33 = de_r_cr(biny, len33, lenbinbit);
					Uint8_Dat z33U8v;
					z33U8v.dat = zCf33.z;
					z33U8v.len = zCf33.lenzbit;
					int lenw3 = sumBin(&z33U8v);
					int lenk3 = zCf33.lenzbit - lenw3;
					DEC cf3k3Dec = de_2n_code_new(biny, lenk3, thdt, typ);
					//Uint32_Dat cf3w3;
					//cf3w3.len = lenw3;
					//cf3w3.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					//for (int i = 0; i < lenw3; i++)
					//{
					//	cf3w3.dat[i] = 1;
					//}
					DEC cf3w3Dec = de_2n_code_new(biny, lenw3, thdt, typ);
					cf23.len = zCf33.lenzbit;
					cf23.dat = separate_inv(zCf33.z, cf3w3Dec.r, cf3k3Dec.r, zCf33.lenzbit, thdt);
					free(zCf33.z);
					free(cf3k3Dec.r);
					free(cf3w3Dec.r);
				}
				//nd3 == 4
				Uint32_Dat cf24;
				cf24.dat = NULL;
				cf24.len = 0;
				if (L3.sym < 4)
				{
					DEC cf22Dec = de_2n_code_new(biny, len34, thdq, typ);
					cf24.dat = cf22Dec.r;
					cf24.len = cf22Dec.lenr;
				}
				else
				{
					//DEC zCf33 = de_r_cr(biny, len33, lenbinbit);
					//cf23.dat = (unsigned int*)calloc(zCf33.lenzbit, sizeof(unsigned int));
					//cf23.len = zCf33.lenzbit;
					//int indexZ1 = 0;
					//for (int i = 0; i < zCf33.lenzbit; i++)
					//{
					//	indexZ1 = i & 7;
					//	if (((zCf33.z[i >> 3] >> (7 - indexZ1)) & 1))
					//	{
					//		cf23.dat[i] = 2;
					//	}
					//	else
					//	{
					//		cf23.dat[i] = 1;
					//	}
					//}
					//free(zCf33.z);
					DEC zCf34 = de_r_cr(biny, len34, lenbinbit);
					Uint8_Dat z34U8v;
					z34U8v.dat = zCf34.z;
					z34U8v.len = zCf34.lenzbit;
					int lenw3 = sumBin(&z34U8v);
					int lenk3 = zCf34.lenzbit - lenw3;
					DEC cf3k4Dec = de_2n_code_new(biny, lenk3, thdt, typ);
					//Uint32_Dat cf3w3;
					//cf3w3.len = lenw3;
					//cf3w3.dat = (unsigned int *)calloc(lenw3, sizeof(unsigned int));
					//for (int i = 0; i < lenw3; i++)
					//{
					//	cf3w3.dat[i] = 1;
					//}
					DEC cf3w4Dec = de_2n_code_new(biny, lenw3, thdt, typ);
					cf24.len = zCf34.lenzbit;
					cf24.dat = separate_inv(zCf34.z, cf3w4Dec.r, cf3k4Dec.r, zCf34.lenzbit, thdt);
					free(zCf34.z);
					free(cf3k4Dec.r);
					free(cf3w4Dec.r);
				}



				cf12.len = z22.lenzbit;
				cf12.dat = separate_inv(z22.z, cf24.dat, cf23.dat, z22.lenzbit, thdq);

				ret.cfk.dat = separate_inv(z1.z, cf12.dat, cf11.dat, z1.lenzbit, thdh);
				ret.cfk.len = z1.lenzbit;

				free(z21.z);
				if (cf21.dat != NULL)
				{
					free(cf21.dat);
				}
				if (cf22.dat != NULL)
				{
					free(cf22.dat);
				}
				if (cf11.dat != NULL)
				{
					free(cf11.dat);
				}
				free(z22.z);
				free(cf23.dat);
				free(cf24.dat);
				free(cf12.dat);
			}
		}
	}
	return ret;
}
