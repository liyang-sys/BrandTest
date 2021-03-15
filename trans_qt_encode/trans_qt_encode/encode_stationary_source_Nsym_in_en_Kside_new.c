#include <stdio.h>
#include "all.h"
#include"encoding.h"

#define uint unsigned int
#define uchar unsigned char

void encode_stationary_source_Nsym_in_en_Kside_new(uint *r, int lenr, float p_or_codebook, int Nsym,int flg)//on_off是否将codebook编进去bin,flg是SFcode的第二个参数
{
	union data {
		unsigned short int a;
		uchar b[4];
	} rem;
	extern int ptr;
	extern uchar *bin;
	int i, x;
	SFC sfc1;

	if (Nsym <= 5)
	{
		//printf("encode_stationary_source_Nsym 1:\n");///
		//ptr++; ///codebook=0
		ptr += 4;//修改的
		x = ptr & 7;
		for (i = 0; i < lenr; i++)
		{
			sfc1 = SFcode(r[i], Nsym);
			rem.a = sfc1.code;
			rem.a = rem.a << (16 - x - sfc1.lb);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr += sfc1.lb; x += sfc1.lb; x &= 7;
		}
	}
	else
	{
		//printf("encode_stationary_source_Nsym 2:\n");///
		//encode_stationary_source_Nsym_sub(r, lenr, p_or_codebook, Nsym, on_off, flg);
		encode_stationary_source_Nsym_sub_in_en_Kside_new(r, lenr, p_or_codebook, Nsym, flg);
	}
}