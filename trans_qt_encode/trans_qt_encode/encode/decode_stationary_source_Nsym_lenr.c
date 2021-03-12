#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char

DEC decode_stationary_source_Nsym_lenr(uchar *bin, int codebook, int lenr0, int Nsym, int lenbinbit)
{
	DEC dec1;
	union data3 {
		uint a;
		uchar b[4];
	} rem;
	int extern ptr;
	int i, x;
	int m1, m2, m, k;
	int a; //存放每次读bit的值
	int mai;
	int hd = 0;
	int hd0;
	int nfrac;
	DES des1;

	dec1.r = (uint *)malloc(lenbinbit * sizeof(uint));

	dec1.lenr = 0;
	x = ptr & 7;
	//m=2类型
	if (codebook << 31 == 0) {
		k = (codebook >> 1) + 1; m = 1 << k;
		if (Nsym > m)
		{
			hd0 = (Nsym - 0.5) / m - 1;
			nfrac = Nsym - hd0 * m;
			while (dec1.lenr < lenr0) {
				mai = 0;
				while ((((bin[ptr >> 3] >> (7 - x)) & 1) == 1) && (hd < hd0))
				{
					mai += m;
					hd++;
					x++; x &= 7; ptr++;
				}
				if (hd < hd0)
				{
					x++; x &= 7; ptr++;
					rem.a = 0;
					rem.b[2] |= bin[ptr >> 3]; rem.b[1] |= bin[(ptr >> 3) + 1]; rem.b[0] |= bin[(ptr >> 3) + 2];
					rem.a = rem.a << (8 + x); rem.a = rem.a >> (32 - k);
					dec1.r[dec1.lenr] = mai + rem.a + 1;
					dec1.lenr++;
					x += k; x &= 7; ptr += k;
				}
				else
				{
					// 处理nfrac==256的极端情况，手动重新解码
					if (nfrac == 256) {
						rem.b[3] = bin[ptr >> 3];
						rem.b[2] = bin[(ptr >> 3) + 1];
						x = ptr & 7;
						des1.sym = ((rem.a << x) >> (32 - 8)) + 1;
						ptr += 8;
					}
					else {
						des1 = deSFcode(bin, nfrac);
					}
					dec1.r[dec1.lenr] = mai + des1.sym;
					dec1.lenr++;
					x = ptr & 7;
				}
				hd = 0;
			}
		}
		else
		{
			while (dec1.lenr < lenr0)
			{
				des1 = deSFcode(bin, Nsym);
				dec1.r[dec1.lenr] = des1.sym;
				dec1.lenr++;
				x = ptr & 7;
			}
		}
	}
	//m=3类型
	else {
		k = codebook >> 1; m1 = 1 << k; m2 = m1 << 1; m = m1 + m2;
		if (Nsym > m1)
		{
			hd0 = (Nsym - 0.5 - m / 3) / m - 1;
			nfrac = Nsym - hd0 * m;
			while (dec1.lenr < lenr0) {
				mai = 0;
				if (hd0 > 0)
				{
					while ((((bin[ptr >> 3] >> (7 - x)) & 1) == 1) && (hd < hd0))
					{
						mai += m;
						hd++;
						x++; x &= 7; ptr++;
					}
					if (hd < hd0)
					{
						x++; x &= 7; ptr++;
						if (k == 0) {
							if (((bin[ptr >> 3] >> (7 - x)) & 1) == 0) {
								dec1.r[dec1.lenr] = mai + 1;
							}
							else {
								x++; x &= 7; ptr++;
								if (((bin[ptr >> 3] >> (7 - x)) & 1) == 0) {
									dec1.r[dec1.lenr] = mai + 2;
								}
								else {
									dec1.r[dec1.lenr] = mai + 3;
								}
							}
							dec1.lenr++;
							x++; x &= 7; ptr++;
						}
						else {
							if (((bin[ptr >> 3] >> (7 - x)) & 1) == 0) {
								x++; x &= 7; ptr++;
								rem.a = 0;
								rem.b[2] |= bin[ptr >> 3]; rem.b[1] |= bin[(ptr >> 3) + 1]; rem.b[0] |= bin[(ptr >> 3) + 2];
								rem.a = rem.a << (8 + x); rem.a = rem.a >> (32 - k);
								dec1.r[dec1.lenr] = mai + rem.a + 1;
							}
							else {
								x++; x &= 7; ptr++;
								if (((bin[ptr >> 3] >> (7 - x)) & 1) == 0) {
									x++; x &= 7; ptr++;
									rem.a = 0;
									rem.b[2] |= bin[ptr >> 3]; rem.b[1] |= bin[(ptr >> 3) + 1]; rem.b[0] |= bin[(ptr >> 3) + 2];
									rem.a = rem.a << (8 + x); rem.a = rem.a >> (32 - k);
									dec1.r[dec1.lenr] = mai + rem.a + 1 + m1;
								}
								else {
									x++; x &= 7; ptr++;
									rem.a = 0;
									rem.b[2] |= bin[ptr >> 3]; rem.b[1] |= bin[(ptr >> 3) + 1]; rem.b[0] |= bin[(ptr >> 3) + 2];
									rem.a = rem.a << (8 + x); rem.a = rem.a >> (32 - k);
									dec1.r[dec1.lenr] = mai + rem.a + 1 + m2;
								}
							}
							dec1.lenr++;
							x += k; x &= 7; ptr += k;
						}
					}
					else
					{
						// 处理nfrac==256的极端情况，手动重新解码
						if (nfrac == 256) {
							rem.b[3] = bin[ptr >> 3];
							rem.b[2] = bin[(ptr >> 3) + 1];
							x = ptr & 7;
							des1.sym = ((rem.a << x) >> (32 - 8)) + 1;
							ptr += 8;
						}
						else {
							des1 = deSFcode(bin, nfrac);
						}
						dec1.r[dec1.lenr] = mai + des1.sym;
						dec1.lenr++;
						x = ptr & 7;
					}
				}
				else
				{
					des1 = deSFcode(bin, Nsym);
					dec1.r[dec1.lenr] = des1.sym;
					dec1.lenr++;
					x = ptr & 7;
				}
				hd = 0;
			}
		}
		else
		{
			while (dec1.lenr < lenr0){
				des1 = deSFcode(bin, Nsym);
				dec1.r[dec1.lenr] = des1.sym;
				dec1.lenr++;
				x = ptr & 7;
			}
		}
	}
	return  dec1;
}
