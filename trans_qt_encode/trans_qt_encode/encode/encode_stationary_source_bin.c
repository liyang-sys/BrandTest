#include <stdio.h>
#include "all.h"
#include<math.h>

#define uint unsigned int
#define uchar unsigned char

void encode_stationary_source_bin(uint *r, int lenr, float pw, int k, int m, int m1, int m2, int cls)
{
	union data {
		uint a;
		uchar b[4];
	} rem;
	extern int ptr;
	extern uchar *bin;
	int i, j, n, x;
	int codebook = 0;

	if (pw >= 1 || pw == 0) {
		codebook = (int)pw;
		if ((codebook % 2) == 0)
		{
			k = (codebook >> 1) + 1; m = 1 << k; cls = 0;
		}
		else
		{
			k = codebook >> 1; m1 = 1 << k; m2 = m1 << 1; m = m1 + m2; cls = 1;
		}
	}

	x = ptr & 7;
	//第一种类型:m=2^k
	if (cls == 0) {
		for (i = 0; i < lenr; i++) {
			rem.a = r[i];
			if (rem.a != 0)
			{
				while (rem.a > m) {
					rem.a -= m;
					bin[ptr >> 3] |= 1 << (7 - x);
					x++; x &= 7; ptr++;
				}
				x++; x &= 7; ptr++;
			}
			//int hd = (int)(floor((rem.a - 0.5) / (m*1.0)));
			//if (hd >= 1)
			//{
			//	for (int i = 0; i < hd; i++)
			//	{
			//			bin[ptr >> 3] |= 1 << (7 - x);
			//			x++; x &= 7; ptr++;
			//	}
			//	x++; x &= 7; ptr++;
			//}
			rem.a = (rem.a - 1) << (24 - x - k);
			bin[ptr >> 3] |= rem.b[2];
			bin[(ptr >> 3) + 1] |= rem.b[1];
			bin[(ptr >> 3) + 2] |= rem.b[0];
			x += k; x &= 7; ptr += k;
		}
	}
	//第二种类型: m=3*2^(k-1)
	else {
		for (i = 0; i < lenr; i++) {
			rem.a = r[i];
			if (rem.a != 0)
			{
				while (rem.a > m) {
					rem.a -= m;
					bin[ptr >> 3] |= 1 << (7 - x);
					x++; x &= 7; ptr++;
				}
				x++; x &= 7; ptr++;
			}

			//int hd = (int)(floor((rem.a - 0.5) / (m*1.0)));
			//if (hd >= 1)
			//{
			//	for (int i = 0; i < hd; i++)
			//	{
			//		bin[ptr >> 3] |= 1 << (7 - x);
			//		x++; x &= 7; ptr++;
			//	}
			//	x++; x &= 7; ptr++;
			//}

			if (rem.a <= m1) {
				x++; x &= 7; ptr++;
			}
			else if (rem.a <= m2) {
				rem.a -= m1;
				bin[ptr >> 3] |= 1 << (7 - x);
				x += 2; x &= 7; ptr += 2;
			}
			else {
				rem.a -= m2;
				bin[ptr >> 3] |= 1 << (7 - x);
				x++; x &= 7; ptr++;
				bin[ptr >> 3] |= 1 << (7 - x);
				x++; x &= 7; ptr++;
			}
			if (k > 0) {
				rem.a = (rem.a - 1) << (24 - x - k);
				bin[ptr >> 3] |= rem.b[2];
				bin[(ptr >> 3) + 1] |= rem.b[1];
				bin[(ptr >> 3) + 2] |= rem.b[0];
				x += k; x &= 7; ptr += k;
			}
		}
	}

	//return codebook;
}
