#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
#include <math.h>


DEC de_2n_code_new(uchar *bin, int len, int maxc, int typ)
{
	extern int ptr;
	DEC c2n;
	c2n.r = calloc(len, sizeof(unsigned int));
	c2n.lenr = len;
	int thd = 0;
	double n = 0;
	int binBytePos = ptr & 7;
	if (typ == 1)
	{
		thd = maxc >> 1;
		n = log(thd) / log(2);
		if (n > 0.5)
		{
			int flg = (bin[ptr >> 3] >> (7 - binBytePos)) & 1;
			ptr++;
			binBytePos &= 7;
			if (flg == 0)
			{
				int thdh = thd >> 1;
				DEC z = de_z0(bin, len);
				int sumZ = 0;
				int zBytePos = 0;
				for (int i = 0; i < z.lenzbit; i++)
				{
					zBytePos = i & 7;
					int temp = ((z.z[i >> 3]) >> (7 - zBytePos)) & 1;
					if (temp == 1)
					{
						sumZ++;
					}
				}

				DEC z1 = de_z0(bin, sumZ);
				int sumZ1 = 0;
				int z1BytePos = 0;
				for (int i = 0; i < z1.lenzbit; i++)
				{
					z1BytePos = i & 7;
					int temp = ((z1.z[i >> 3]) >> (7 - z1BytePos)) & 1;
					if (temp == 1)
					{
						sumZ1++;
					}
				}

				DEC z2 = de_z0(bin, sumZ1);
				Uint8_Dat z2Uint8;
				z2Uint8.dat = z2.z;
				z2Uint8.len = z2.lenzbit;
				Uint32_Dat ckUint32, clkUint32, c2Uint32;
				if (n > 1.5)
				{
					DEC ck = de_2n_code_sub1(bin, n -1, z.lenzbit - sumZ);
					ckUint32.dat = ck.r;
					ckUint32.len = ck.lenr;

					DEC clk = de_2n_code_sub1(bin, n - 1, z1.lenzbit - sumZ1);
					for (int i = 0; i < clk.lenr; i++)
					{
						clk.r[i] += thdh;
					}
					clkUint32.dat = clk.r;
					clkUint32.len = clk.lenr;

					DEC c2 = de_2n_code_sub(bin,n-1,&z2Uint8);
					for (int i = 0; i < c2.lenr; i++)
					{
						c2.r[i] += thd;
					}
					c2Uint32.dat = c2.r;
					c2Uint32.len = c2.lenr;
				}
				else
				{
					ckUint32.len = z.lenzbit - sumZ;
					ckUint32.dat = calloc(ckUint32.len, sizeof(unsigned int));
					for (int i = 0; i < ckUint32.len; i++)
					{
						ckUint32.dat[i] = 1;
					}

					clkUint32.len = z1.lenzbit - sumZ1;
					clkUint32.dat = calloc(clkUint32.len, sizeof(unsigned int));
					for (int i = 0; i < clkUint32.len; i++)
					{
						clkUint32.dat[i] = 2;
					}

					c2Uint32.len = z2.lenzbit;
					c2Uint32.dat = calloc(c2Uint32.len, sizeof(unsigned int));
					int z2BytePos = 0;
					int temp = 0;
					for (int i = 0; i < z2.lenzbit; i++)
					{
						z2BytePos = i & 7;
						temp = ((z2.z[i >> 3]) >> (7 - z2BytePos)) & 1;
						c2Uint32.dat[i] = temp + 3;
					}
				}
				unsigned int* c1 = separate_inv(z1.z, c2Uint32.dat, clkUint32.dat, z1.lenzbit, 0);
				c2n.r = separate_inv(z.z, c1, ckUint32.dat, z.lenzbit, 0);
				c2n.lenr = z.lenzbit;

				if (z.z != NULL)
				{
					free(z.z);
				}
				if (z1.z != NULL)
				{
					free(z1.z);
				}
				if (z2.z != NULL)
				{
					free(z2.z);
				}
				if (ckUint32.dat!= NULL)
				{
					free(ckUint32.dat);
				}
				if (clkUint32.dat != NULL)
				{
					free(clkUint32.dat);
				}
				if (c2Uint32.dat != NULL)
				{
					free(c2Uint32.dat);
				}
				if (c1 != NULL)
				{
					free(c1);
				}
			}
			else
			{
				DEC ck,cw;
				DEC z = de_z0(bin, len);
				//printf("de_z0ºóptr = %d\n", ptr);
				int lenw = 0;
				int zBytePos = 0;
				for (int i = 0; i < z.lenzbit; i++)
				{
					zBytePos = i & 7;
					int temp = ((z.z[i >> 3]) >> (7 - zBytePos)) & 1;
					if (temp == 1)
					{
						lenw++;
					}
				}
				int lenk = len - lenw;
				DEC zk = de_z0(bin, lenk);
				DEC zw = de_z0(bin, lenw);
				if (n > 1.5)
				{
					Uint8_Dat zkUint8;
					zkUint8.dat = zk.z;
					zkUint8.len = zk.lenzbit;
					ck = de_2n_code_sub(bin, n - 1, &zkUint8);

					Uint8_Dat zwUint8;
					zwUint8.dat = zw.z;
					zwUint8.len = zw.lenzbit;
					cw = de_2n_code_sub(bin, n - 1, &zwUint8);

				}
				else
				{
					ck.lenr = zk.lenzbit;
					ck.r = calloc(ck.lenr, sizeof(unsigned int));
					int zkBytePos = 0;
					int temp = 0;
					for (int i = 0; i < zk.lenzbit; i++)
					{
						zkBytePos = i & 7;
						temp = ((zk.z[i >> 3]) >> (7 - zkBytePos)) & 1;
						ck.r[i] = temp + 1;
					}

					cw.lenr = zw.lenzbit;
					cw.r = calloc(cw.lenr, sizeof(unsigned int));
					int zwBytePos = 0;
					for (int i = 0; i < zw.lenzbit; i++)
					{
						zwBytePos = i & 7;
						temp = ((zw.z[i >> 3]) >> (7 - zwBytePos)) & 1;
						cw.r[i] = temp + 1;
					}
				}
				c2n.r = separate_inv(z.z, cw.r, ck.r, z.lenzbit, thd);
				c2n.lenr = z.lenzbit;
				if (ck.r != NULL)
				{
					free(ck.r);
				}
				if (cw.r != NULL)
				{
					free(cw.r);
				}
				if (zk.z != NULL)
				{
					free(zk.z);
				}
				if (zw.z != NULL)
				{
					free(zw.z);
				}
				if (z.z != NULL)
				{
					free(z.z);
				}
			}
		}
		else
		{
			DEC z = de_z0(bin, len);
			int zBytePos = 0;
			for (int i = 0; i < z.lenzbit; i++)
			{
				zBytePos = i & 7;
				int temp = ((z.z[i >> 3]) >> (7 - zBytePos)) & 1;
				c2n.r[i] = temp + 1;
			}
			if (z.z != NULL)
			{
				free(z.z);
			}
		}
	}
	else
	{
		thd = maxc /3;
		n = log(thd) / log(2);
		DEC ck, cw;
		DEC z = de_z0(bin, len);
		int lenw = 0;
		int zBytePos = 0;
		for (int i = 0; i < z.lenzbit; i++)
		{
			zBytePos = i & 7;
			int temp = ((z.z[i >> 3]) >> (7 - zBytePos)) & 1;
			if (temp == 1)
			{
				lenw++;
			}
		}
		int lenk = len - lenw;
		if (n > 0.5)
		{
			DEC zk = de_z0(bin, lenk);
			DEC zw = de_z0(bin, lenw);
			if (n > 1.5)
			{
				Uint8_Dat zkUint8;
				zkUint8.dat = zk.z;
				zkUint8.len = zk.lenzbit;
				ck = de_2n_code_sub(bin, n - 1, &zkUint8);

				Uint8_Dat zwUint8;
				zwUint8.dat = zw.z;
				zwUint8.len = zw.lenzbit;
				cw = de_2n_code_sub(bin, n, &zwUint8);
				//for (int i = 0; i < cw.lenr; i++)
				//{
				//	cw.r[i] += thd;
				//}
			}
			else
			{
				ck.lenr = zk.lenzbit;
				ck.r = calloc(ck.lenr, sizeof(unsigned int));
				int zkBytePos = 0;
				int temp = 0;
				for (int i = 0; i < zk.lenzbit; i++)
				{
					zkBytePos = i & 7;
					temp = ((zk.z[i >> 3]) >> (7 - zkBytePos)) & 1;
					ck.r[i] = temp + 1;
				}

				Uint8_Dat zwUint8;
				zwUint8.dat = zw.z;
				zwUint8.len = zw.lenzbit;
				cw = de_2n_code_sub(bin, n, &zwUint8);
			}
			c2n.r = separate_inv(z.z, cw.r, ck.r, z.lenzbit, thd);
			c2n.lenr = z.lenzbit;
			if (zk.z != NULL)
			{
				free(zk.z);
			}
			if (zw.z != NULL)
			{
				free(zw.z);
			}
		}
		else
		{
			DEC zw = de_z0(bin, lenw);
			cw.lenr = zw.lenzbit;
			cw.r = calloc(cw.lenr, sizeof(unsigned int));
			int zwBytePos = 0;
			for (int i = 0; i < zw.lenzbit; i++)
			{
				zwBytePos = i & 7;
				cw.r[i] = ((zw.z[i >> 3]) >> (7 - zwBytePos)) & 1;
			}

			ck.lenr = lenk;
			ck.r = calloc(ck.lenr, sizeof(unsigned int));
			for (int i = 0; i < lenk; i++)
			{
				ck.r[i] = 1;
			}
			c2n.r = separate_inv(z.z, cw.r, ck.r, z.lenzbit, 2);
			c2n.lenr = z.lenzbit;
			if (zw.z != NULL)
			{
				free(zw.z);
			}
		}
		if (ck.r != NULL)
		{
			free(ck.r);
		}
		if (cw.r != NULL)
		{
			free(cw.r);
		}
		if (z.z != NULL)
		{
			free(z.z);
		}
	}
	return c2n;
}