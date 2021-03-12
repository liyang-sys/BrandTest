#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
#include <math.h>


void en_2n_code_new(Uint32_Dat* c2n, Uint32_Dat* nc,int typ,int maxc)
{
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
    extern int ptr;
    extern uchar *bin;
    int thd = 0;
    double n = 0;
    int thdh = 0;
    if(typ == 1)
    {
        thd = maxc >> 1;
        n = log(thd)/log(2);
        if(n > 0.5)
        {
            thdh = thd >> 1;
            int sumNc = sum(nc->dat,nc->len);
            double p1 = sum(nc->dat,thdh)/(double)sumNc;
            double p2 = sum(nc->dat,thd)/(double)sumNc;
            if(fabs(p1 - 0.5) < fabs(p2 - 0.5))
            {
				int binBytePos = ptr & 7;
				ptr++; binBytePos++;
				binBytePos &= 7;//编码0
				SEP sepC2n = separate0(c2n->dat, thdh, c2n->len, 1);//c2n->dat是无符号的int
				Uint8_Dat z;
				z.dat = sepC2n.sep;
				z.len = sepC2n.lensepbit;//有效位
				en_z0(&z, sepC2n.lrw); //编码binz

				SEP sepC1 = separate0(sepC2n.rw, thdh, sepC2n.lrw, 1);
				Uint8_Dat z1;
				z1.dat = sepC1.sep;
				z1.len = sepC1.lensepbit;//有效位
				en_z0(&z1, sepC1.lrw);//编码binz1

				SEP sepC2 = separate0(sepC1.rw, thdh, sepC1.lrw, 1);
				Uint8_Dat z2;
				z2.dat = sepC2.sep;
				z2.len = sepC2.lensepbit;//有效位
				en_z0(&z2, sepC2.lrw);//编码binz2

				Uint32_Dat ck;
				ck.dat = c2n->dat;
				ck.len = sepC2n.lrk;
				en_2n_code_sub(&ck, NULL, n - 1);//编码bink

				Uint32_Dat clk;
				clk.dat = sepC2n.rw;
				clk.len = sepC1.lrk;
				en_2n_code_sub(&clk, NULL, n - 1);//编码bink1

				Uint32_Dat c3;
				c3.dat = sepC2.rw;
				c3.len = sepC2.lrw;

				Uint32_Dat c2k;
				c2k.dat = sepC1.rw;
				c2k.len = sepC2.lrk;
				en_2n_code_sub(&c3, &c2k, n - 1);//编码bin2

				if (sepC2n.rw != NULL)
				{
					free(sepC2n.rw);
				}
				if (sepC2n.sep != NULL)
				{
					free(sepC2n.sep);
				}
				if (sepC1.rw != NULL)
				{
					free(sepC1.rw);
				}
				if (sepC1.sep != NULL)
				{
					free(sepC1.sep);
				}
				if (sepC2.rw != NULL)
				{
					free(sepC2.rw);
				}
				if (sepC2.sep != NULL)
				{
					free(sepC2.sep);
				}

            }
            else
            {
				int binBytePos = ptr & 7; //编码1
				rem.a = 1;
				rem.a = rem.a << (15 - binBytePos);
				bin[ptr >> 3] |= rem.b[1];
				bin[(ptr >> 3) + 1] |= rem.b[0];
				ptr ++; binBytePos++;
				binBytePos &= 7;

                SEP sepC2n = separate0(c2n->dat,thd,c2n->len,1);//c2n->dat是无符号的int
                Uint8_Dat z ;
                z.dat =  sepC2n.sep;
                z.len = sepC2n.lensepbit;//有效位
                en_z0(&z,sepC2n.lrw); //编码binz

				SEP sepCk = separate0(c2n->dat, thdh, sepC2n.lrk, 1);
				Uint8_Dat zk;
				zk.dat = sepCk.sep;
				zk.len = sepCk.lensepbit;//有效位




				en_z0(&zk, sepCk.lrw);//编码binzk

				//FILE*fp = fopen("binTC.txt", "wb");
				//fwrite(bin, sizeof(unsigned char), (ptr / 8) + 1, fp);
				//fclose(fp);

				SEP sepCw = separate0(sepC2n.rw, thdh, sepC2n.lrw, 1);
				Uint8_Dat zw;
				zw.dat = sepCw.sep;
				zw.len = sepCw.lensepbit;//有效位
				en_z0(&zw, sepCw.lrw);//编码binzw

				Uint32_Dat ckw;
				ckw.dat = sepCk.rw;
				ckw.len = sepCk.lrw;

				Uint32_Dat ckk;
				ckk.dat = c2n->dat;
				ckk.len = sepCk.lrk;
				en_2n_code_sub(&ckw, &ckk, n - 1);//编码bin2k

				Uint32_Dat cww;
				cww.dat = sepCw.rw;
				cww.len = sepCw.lrw;

				Uint32_Dat cwk;
				cwk.dat = sepC2n.rw;
				cwk.len = sepCw.lrk;
				en_2n_code_sub(&cww, &cwk, n - 1);//编码bin2w
				//printf("\n编码结束的ptr = %d\n", ptr);
    //            int temp = 0;
				//释放内存
				if (sepC2n.rw != NULL)
				{
					free(sepC2n.rw);
				}
				if (sepC2n.sep != NULL)
				{
					free(sepC2n.sep);
				}
				if (sepCk.rw != NULL)
				{
					free(sepCk.rw);
				}
				if (sepCk.sep != NULL)
				{
					free(sepCk.sep);
				}
				if (sepCw.rw != NULL)
				{
					free(sepCw.rw);
				}
				if (sepCw.sep != NULL)
				{
					free(sepCw.sep);
				}
            }
        }
        else
        {
			SEP sepC2n = separate0(c2n->dat, thd, c2n->len, 1);//c2n->dat是无符号的int
			Uint8_Dat z;
			z.dat = sepC2n.sep;
			z.len = sepC2n.lensepbit;//有效位
			en_z0(&z, sepC2n.lrw); //编码binz
			if (sepC2n.rw != NULL)
			{
				free(sepC2n.rw);
			}
			if (sepC2n.sep != NULL)
			{
				free(sepC2n.sep);
			}
			//printf("\n en_z0编码 z 结束的ptr = %d\n", ptr);
        }
    }
    else
    {
        thd = maxc/3;
        n = log(thd)/log(2);
        SEP sepC2n = separate0(c2n->dat,thd,c2n->len,1);//c2n->dat是无符号的int
        Uint8_Dat z ;
        z.dat =  sepC2n.sep;
        z.len = sepC2n.lensepbit;//有效位
        en_z0(&z,sepC2n.lrw);
		//printf("\n en_z0编码 z 结束的ptr = %d\n", ptr);
		if (n > 0.5)
		{
			thdh = thd >> 1;

			SEP sepCk = separate0(c2n->dat, thdh, sepC2n.lrk, 1);
			Uint8_Dat zk;
			zk.dat = sepCk.sep;
			zk.len = sepCk.lensepbit;//有效位
			en_z0(&zk, sepCk.lrw);//编码binzk
			//printf("\n en_z0编码 zk 结束的ptr = %d\n", ptr);

			SEP sepCw = separate0(sepC2n.rw, thd, sepC2n.lrw, 1);
			Uint8_Dat zw;
			zw.dat = sepCw.sep;
			zw.len = sepCw.lensepbit;//有效位
			en_z0(&zw, sepCw.lrw);//编码binzw
			//printf("\n en_z0编码 zw结束的ptr = %d\n", ptr);

			Uint8_Dat ckw;
			ckw.dat = sepCk.rw;
			ckw.len = sepCk.lrw;

			Uint8_Dat ckk;
			ckk.dat = c2n->dat;
			ckk.len = sepCk.lrk;
			en_2n_code_sub(&ckw, &ckk, n - 1);//编码bin2k

			Uint8_Dat cww;
			cww.dat = sepCw.rw;
			cww.len = sepCw.lrw;

			Uint8_Dat cwk;
			cwk.dat = sepC2n.rw;
			cwk.len = sepCw.lrk;
			en_2n_code_sub(&cww, &cwk, n);//编码bin2w
			//printf("\n编码结束的ptr = %d\n", ptr);
			//释放内存
			if (sepCk.rw != NULL)
			{
				free(sepCk.rw);
			}
			if (sepCk.sep != NULL)
			{
				free(sepCk.sep);
			}
			if (sepCw.rw != NULL)
			{
				free(sepCw.rw);
			}
			if (sepCw.sep != NULL)
			{
				free(sepCw.sep);
			}
		}
		else
		{
			int index = 0;
			int indexZ = 0;
			int temp = 0;
			Uint8_Dat cw;
			cw.dat = calloc((sepC2n.lrw / 8 + 1), sizeof(unsigned char));
			cw.len = sepC2n.lrw;
			int wuns = 0;
			for (int i = 0; i < cw.len; i++)
			{
				temp = sepC2n.rw[i] - 1;
				if (temp == 1)
				{
					rem.a = 1;
					rem.a = rem.a << (15 - index);
					cw.dat[indexZ >> 3] |= rem.b[1];
					cw.dat[(indexZ >> 3) + 1] |= rem.b[0];
					wuns++;
				}
				index++;
				indexZ++;
				index &= 7;
			}
			en_z0(&cw, wuns);
			if (cw.dat != NULL)
			{
				free(cw.dat);
			}
			//printf("\n en_z0编码 z 结束的ptr = %d\n", ptr);
		}
		if (sepC2n.rw != NULL)
		{
			free(sepC2n.rw);
		}
		if (sepC2n.sep != NULL)
		{
			free(sepC2n.sep);
		}
		
    }
}
