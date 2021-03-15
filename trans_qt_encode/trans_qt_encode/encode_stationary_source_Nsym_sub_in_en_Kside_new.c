#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char


void encode_stationary_source_Nsym_sub_in_en_Kside_new(uint *r, int lenr, float pw, int Nsym, int flg)
{
	extern int ptr;
	extern uchar *bin;
	union data {
		uint a;
		uchar b[4];
	} rem;
	int codebook;
	int num = 1;
	int i, j, x, m1, m2, m, k, cls; //m1表示m/3,m2表示2*m/3
	SFC se, sfcode;
	int n;
	int hd = 0;
	int hd0;
	int nfrac;
	int t, t0, t1;
	//e.len_bit = 0;



	if (pw < 1 && pw != 0) {
		if (pw > 0.242141716744801)
			codebook = 0, k = 1, m = 2, cls = 0; //m=2
		else if (pw > 0.179664643992362)
			codebook = 1, k = 0, m = 3, m1 = 1, m2 = 2, cls = 1; //m=3
		else if (pw > 0.129449436703876)
			codebook = 2, k = 2, m = 4, cls = 0; //m=4
		else if (pw > 0.094276335736093)
			codebook = 3, k = 1, m = 6, m1 = 2, m2 = 4, cls = 1; //m=6

		else if (pw > 0.066967008463193)
			codebook = 4, k = 3, m = 8, cls = 0; //m=8
		else if (pw > 0.048304846989380)
			codebook = 5, k = 2, m = 12, m1 = 4, m2 = 8, cls = 1; //m=12
		else if (pw > 0.034063671075154)
			codebook = 6, k = 4, m = 16, cls = 0; //m=16
		else if (pw > 0.024451357947427)
			codebook = 7, k = 3, m = 24, m1 = 8, m2 = 16, cls = 1; //m=24

		else if (pw > 0.017179401454749)
			codebook = 8, k = 5, m = 32, cls = 0; //m=32
		else if (pw > 0.012301340462298)
			codebook = 9, k = 4, m = 48, m1 = 16, m2 = 32, cls = 1; //m=48
		else if (pw > 0.008626912537338)
			codebook = 10, k = 6, m = 64, cls = 0; //m=64
		else if (pw > 0.006169702847764)
			codebook = 11, k = 5, m = 96, m1 = 32, m2 = 64, cls = 1; //m=96
		else if (pw > 0.004322799566716)
			codebook = 12, k = 7, m = 128, cls = 0; //m=128
		else if (pw > 0.003089624313080)
			codebook = 13, k = 6, m = 192, m1 = 64, m2 = 128, cls = 1; //m=192
		else if (pw > 0.002163740670202)
			codebook = 14, k = 8, m = 256, cls = 0; //m=256
		else if (pw > 0.001546007225711)
			codebook = 15, k = 7, m = 384, m1 = 128, m2 = 256, cls = 1; //m=384

		else
			codebook = 16, k = 9, m = 512, cls = 0; //m=131072
	}
	else
	{
		codebook = (int)pw;
		//if ((codebook % 2) == 0)
		//{
		//	k = (codebook >> 1) + 1; m = 1 << k; cls = 0;
		//}
		//else
		//{
		//	k = codebook >> 1; m1 = 1 << k; m2 = m1 << 1; m = m1 + m2; cls = 1;
		//}
	}

	//printf("codebook = %d cls = %d\n", codebook, cls);

	//if (!on_off)
	//{
	//	for (i = 3; i >= 0; i--)
	//	{
	//		t0 = (codebook >> i) & 1;  ///codebook的值直接赋给bin
	//		t1 = ptr & 7;
	//		bin[ptr >> 3] |= (t0 << (7 - t1));
	//		ptr++;
	//	}
	//}
	//else
	//{
	//	sfcode = SFcode(codebook + 1, flg);
	//	//printf("sfcode.lb = %d\n", sfcode.lb);
	//	for (i = sfcode.lb - 1; i >= 0; i--)
	//	{
	//		t0 = (sfcode.code >> i) & 1;  ///SFcode(codebook+1,16)赋给bin
	//		t1 = ptr & 7;
	//		bin[ptr >> 3] |= (t0 << (7 - t1));
	//		ptr++;
	//	}
	//}
	uint*rTemp = calloc(1, sizeof(uint));
	rTemp[0] = codebook + 1;
	encode_stationary_source_bin(rTemp, 1, flg, 0, 0, 0, 0, 0);

	int ptr1 = ptr;

	x = ptr & 7;

	//第一种类型:m=2^k
	if (!cls) {
		if (Nsym > m)
		{
			hd0 = (Nsym - 0.5) / m - 1;
			nfrac = Nsym - hd0 * m;

			for (j = 0; j < lenr; j++)
			{

				hd = (r[j] - 0.5) / m;
				// 修改
				if (hd > hd0)    hd = hd0;
				//printf("(%d):%d\n", j+1,hd);
				if (hd < hd0)
				{
					//printf("1.1 - (%d): lenb = %d\n", j+1 ,ptr-ptr_org);
					rem.a = r[j];
					while (rem.a > m) {
						rem.a -= m;
						bin[ptr >> 3] |= 1 << (7 - x);
						x++; x &= 7; ptr++;
					}
					x++; x &= 7; ptr++;
					rem.a = (rem.a - 1) << (24 - x - k);
					bin[ptr >> 3] |= rem.b[2];
					bin[(ptr >> 3) + 1] |= rem.b[1];
					bin[(ptr >> 3) + 2] |= rem.b[0];
					x += k; x &= 7; ptr += k;
					//printf("1.2 - (%d): lenb = %d\n", j+1 ,ptr-ptr_org);
				}
				else
				{
					//printf("2.1 - (%d): lenb = %d  r[j]-hd0*m = %d nfrac = %d \n", j+1 ,ptr-ptr_org, r[j] - hd0 * m, nfrac);
					hd = hd0;
					while (hd > 0) {
						hd--;
						bin[ptr >> 3] |= 1 << (7 - x);
						x++; x &= 7; ptr++;
					}
					se = SFcode((r[j] - hd0 * m), nfrac);
					//if (nfrac == 256) {
					//	se.code = r[j] - 1;
					//	se.lb = 8;
					//}
					//printf("len = %d code = 0x%x\n", se.lb, se.code);
					if (se.lb <= 8)
					{
						rem.a = se.code;
						rem.a = rem.a << (16 - x - se.lb);
						bin[ptr >> 3] |= rem.b[1];
						bin[(ptr >> 3) + 1] |= rem.b[0];
					}
					else if (se.lb <= 16 && se.lb > 8)
					{
						//printf("SFcode超过了8位 in encode_stationary_source_Nsym_sub\n");
						rem.a = se.code;
						rem.a = rem.a << (24 - x - se.lb);
						bin[ptr >> 3] |= rem.b[2];
						bin[(ptr >> 3) + 1] |= rem.b[1];
						bin[(ptr >> 3) + 2] |= rem.b[0];
					}
					else if (se.lb <= 24 && se.lb > 16)
					{
						//printf("SFcode超过了24位 in encode_stationary_source_Nsym_sub\n");
						rem.a = se.code;
						rem.a = rem.a << (32 - x - se.lb);
						bin[ptr >> 3] |= rem.b[3];
						bin[(ptr >> 3) + 1] |= rem.b[2];
						bin[(ptr >> 3) + 2] |= rem.b[1];
						bin[(ptr >> 3) + 3] |= rem.b[0];
					}
					else
					{
						printf("SFcode超过了32位 in encode_stationary_source_Nsym_sub\n");
					}

					//rem.a = se.code;
					//rem.a = rem.a << (16 - x - se.lb);
					//bin[ptr >> 3] |= rem.b[1];
					//bin[(ptr >> 3) + 1] |= rem.b[0];
					ptr += se.lb; x += se.lb; x &= 7;
					//printf("2.2 - (%d): lenb = %d lenc = %d\n", j+1 ,ptr-ptr_org, se.lb);
				}
			}
		}
		else
		{
			for (j = 0; j < lenr; j++)
			{
				se = SFcode(r[j], Nsym);
				//if (Nsym == 256) {
				//	se.code = r[j] - 1;
				//	se.lb = 8;
				//}
				//rem.a = se.code;
				//rem.a = rem.a << (16 - x - se.lb);
				//bin[ptr >> 3] |= rem.b[1];
				//bin[(ptr >> 3) + 1] |= rem.b[0];
				if (se.lb <= 8)
				{
					rem.a = se.code;
					rem.a = rem.a << (16 - x - se.lb);
					bin[ptr >> 3] |= rem.b[1];
					bin[(ptr >> 3) + 1] |= rem.b[0];
				}
				else if (se.lb <= 16 && se.lb > 8)
				{
					//printf("SFcode超过了8位 in encode_stationary_source_Nsym_sub\n");
					rem.a = se.code;
					rem.a = rem.a << (24 - x - se.lb);
					bin[ptr >> 3] |= rem.b[2];
					bin[(ptr >> 3) + 1] |= rem.b[1];
					bin[(ptr >> 3) + 2] |= rem.b[0];
				}
				else if (se.lb <= 24 && se.lb > 16)
				{
					//printf("SFcode超过了24位 in encode_stationary_source_Nsym_sub\n");
					rem.a = se.code;
					rem.a = rem.a << (32 - x - se.lb);
					bin[ptr >> 3] |= rem.b[3];
					bin[(ptr >> 3) + 1] |= rem.b[2];
					bin[(ptr >> 3) + 2] |= rem.b[1];
					bin[(ptr >> 3) + 3] |= rem.b[0];
				}
				else
				{
					printf("SFcode超过了32位 in encode_stationary_source_Nsym_sub\n");
				}
				ptr += se.lb; x += se.lb; x &= 7;
			}
		}
	}
	else
	{
		/* hd0=floor((Nsym-0.5-n)/(3*n)) -1; */
		hd0 = (Nsym - 0.5 - m / 3) / m - 1;
		/* nfrac = Nsym-hd0*3*n; */
		nfrac = Nsym - hd0 * m;
		if (hd0 > 0)
		{
			/* for ia=1:length(r) */
			for (j = 0; j < lenr; j++)
			{
				/* hd=floor((r(ia)-0.5)/(3*n)); */
				hd = (r[j] - 0.5) / m;
				/* if hd>hd0, hd=hd0; end */
				if (hd > hd0) hd = hd0;
				if (hd < hd0)
				{
					rem.a = r[j];
					while (rem.a > m) {
						rem.a -= m;
						bin[ptr >> 3] |= 1 << (7 - x);
						x++; x &= 7; ptr++;
					}
					x++; x &= 7; ptr++;
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
					//printf("(1)bin len = %d\n", ptr-ptr1);
				}
				else
				{
					hd = hd0;
					while (hd > 0) {
						hd--;
						bin[ptr >> 3] |= 1 << (7 - x);
						x++; x &= 7; ptr++;
					}
					se = SFcode((r[j] - hd0 * m), nfrac);
					//if (nfrac == 256) {
					//	se.code = r[j] - 1;
					//	se.lb = 8;
					//}
					//rem.a = se.code;
					//rem.a = rem.a << (16 - x - se.lb);
					//bin[ptr >> 3] |= rem.b[1];
					//bin[(ptr >> 3) + 1] |= rem.b[0];
					if (se.lb <= 8)
					{
						rem.a = se.code;
						rem.a = rem.a << (16 - x - se.lb);
						bin[ptr >> 3] |= rem.b[1];
						bin[(ptr >> 3) + 1] |= rem.b[0];
					}
					else if (se.lb <= 16 && se.lb > 8)
					{
						//printf("SFcode超过了8位 in encode_stationary_source_Nsym_sub\n");
						rem.a = se.code;
						rem.a = rem.a << (24 - x - se.lb);
						bin[ptr >> 3] |= rem.b[2];
						bin[(ptr >> 3) + 1] |= rem.b[1];
						bin[(ptr >> 3) + 2] |= rem.b[0];
					}
					else if (se.lb <= 24 && se.lb > 16)
					{
						//printf("SFcode超过了24位 in encode_stationary_source_Nsym_sub\n");
						rem.a = se.code;
						rem.a = rem.a << (32 - x - se.lb);
						bin[ptr >> 3] |= rem.b[3];
						bin[(ptr >> 3) + 1] |= rem.b[2];
						bin[(ptr >> 3) + 2] |= rem.b[1];
						bin[(ptr >> 3) + 3] |= rem.b[0];
					}
					else
					{
						printf("SFcode超过了32位 in encode_stationary_source_Nsym_sub\n");
					}
					ptr += se.lb; x += se.lb; x &= 7;
					//printf("   (2)code = %d len = %d nfrac=%d\n", (r[j] - hd0 * m), se.lb, nfrac);
					//printf("(2)bin len = %d\n", ptr-ptr1);
				}
			}

		}
		else
		{
			for (j = 0; j < lenr; j++)
			{
				se = SFcode(r[j], Nsym);
				//if (Nsym == 256) {
				//	se.code = r[j] - 1;
				//	se.lb = 8;
				//}
				//rem.a = se.code;
				//rem.a = rem.a << (16 - x - se.lb);
				//bin[ptr >> 3] |= rem.b[1];
				//bin[(ptr >> 3) + 1] |= rem.b[0];
				if (se.lb <= 8)
				{
					rem.a = se.code;
					rem.a = rem.a << (16 - x - se.lb);
					bin[ptr >> 3] |= rem.b[1];
					bin[(ptr >> 3) + 1] |= rem.b[0];
				}
				else if (se.lb <= 16 && se.lb > 8)
				{
					//printf("SFcode超过了8位 in encode_stationary_source_Nsym_sub\n");
					rem.a = se.code;
					rem.a = rem.a << (24 - x - se.lb);
					bin[ptr >> 3] |= rem.b[2];
					bin[(ptr >> 3) + 1] |= rem.b[1];
					bin[(ptr >> 3) + 2] |= rem.b[0];
				}
				else if (se.lb <= 24 && se.lb > 16)
				{
					//printf("SFcode超过了24位 in encode_stationary_source_Nsym_sub\n");
					rem.a = se.code;
					rem.a = rem.a << (32 - x - se.lb);
					bin[ptr >> 3] |= rem.b[3];
					bin[(ptr >> 3) + 1] |= rem.b[2];
					bin[(ptr >> 3) + 2] |= rem.b[1];
					bin[(ptr >> 3) + 3] |= rem.b[0];
				}
				else
				{
					printf("SFcode超过了32位 in encode_stationary_source_Nsym_sub\n");
				}
				ptr += se.lb; x += se.lb; x &= 7;
			}
		}
	}


	//while(1);

}

#if 0
if (hd0 > 0)
{
	for (j = 0; j < lenr; j++)
	{

		hd = (r[j] - 0.5) / m;
		if (hd < hd0)
		{
			rem.a = r[j];
			while (rem.a > m) {
				rem.a -= m;
				bin[ptr >> 3] |= 1 << (7 - x);
				x++; x &= 7; ptr++;
			}
			x++; x &= 7; ptr++;
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
		else
		{
			hd = hd0;
			while (hd > 0) {
				hd--;
				bin[ptr >> 3] |= 1 << (7 - x);
				x++; x &= 7; ptr++;
			}
			se = SFcode((r[j] - hd0 * m), nfrac);
			if (nfrac == 256) {
				se.code = r[j] - 1;
				se.lb = 8;
			}
			rem.a = se.code;
			rem.a = rem.a << (16 - x - se.lb);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr += se.lb; x += se.lb; x &= 7;
		}
	}
}
else
{
	for (j = 0; j < lenr; j++)
	{
		se = SFcode(r[j], Nsym);
		if (Nsym == 256) {
			se.code = r[j] - 1;
			se.lb = 8;
		}
		rem.a = se.code;
		rem.a = rem.a << (16 - x - se.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += se.lb; x += se.lb; x &= 7;
	}
}
#endif
