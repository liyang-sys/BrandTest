#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

void en_r_cr( Uint32_Dat* r, int lB,Uint32_Dat* nr,float sumcr, float sumcrc, float pr,Uint8_Dat* z)
{
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
	extern int ptr;
	extern uchar *bin;
	int binBytePos = ptr & 7;
    //int sumz = 0; //◊¢“‚£∫≤‚ ‘µƒ ±∫Ú◊¢ Õ¡À
    //for(int i = 0; i < z->len;i++)
    //{
    //    sumz += z->dat[i];
    //}
	int sumz = sumBin(z);
	if ((sumcr > 90 || sumcrc >= 36) && sumz > 600)
	{
		rem.a = 3;//±‡11
		rem.a = rem.a << (14 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += 2; binBytePos += 2;
		
		binBytePos &= 7;
		rem.a = lB;//±‡lB
		rem.a = rem.a << (15 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr++; binBytePos++;

		en_zone_new(r, nr);
	}
	else if ((sumcr > 39 || sumcrc >= 16) && sumz > 200)
	{
		//ptr = 0;//test”√
		//bin = calloc(10240, sizeof(uchar));
		//binBytePos = ptr & 7;

		rem.a = 2;//±‡10
		rem.a = rem.a << (14 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += 2; binBytePos += 2;
		binBytePos &= 7;

		rem.a = lB;//±‡lB
		rem.a = rem.a << (15 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr++; binBytePos++;
		binBytePos &= 7;

		en_zone_sub2(r);
	}
	else if ((sumcr > 16 || sumcrc >= 7) && (pr<0.43 | pr>0.57) && (sumz > 60))
	{
		rem.a = 1;//±‡01
		rem.a = rem.a << (14 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += 2; binBytePos += 2;
		binBytePos &= 7;

		rem.a = lB;//±‡lB
		rem.a = rem.a << (15 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr++; binBytePos++;
		binBytePos &= 7;
		int *n = NULL;
		en_r1(r->dat, 0.8125, r->len, n,0);
	}
	else
	{
		ptr += 2;
		en_z0_r0(z, r, lB, sumBin(z), pr);
	}
}
