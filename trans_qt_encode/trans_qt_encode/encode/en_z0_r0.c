#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

void en_z0_r0(Uint8_Dat* z, Uint32_Dat* runs,int lastBit,int wuns,float p)//wuns是z的和
{
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
	extern int ptr;
	extern uchar *bin;
	double f0 = 0.4219, f1 = 0.3750, f2 = 0.3125;
	int binBytePos = ptr & 7;
	if (wuns >= (z->len) *f0)
	{
		SFC sfc = SFcode(1, 36);
		rem.a = sfc.code;
		rem.a = rem.a << (16 - binBytePos - sfc.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc.lb; binBytePos += sfc.lb; binBytePos &= 7;

		int t1, t0, t;
		for (int i = 0; i < z->len; i++) {   //copy彭前师兄的
			t1 = i & 7;
			t0 = ptr & 7;
			t = ((z->dat[i >> 3]) >> (7 - t1)) & 1;
			bin[ptr >> 3] |= (t << (7 - t0));
			ptr++;
		}
	}
	else if (wuns > (z->len) *f1)
	{
		SFC sfc = SFcode(2, 36);
		rem.a = sfc.code;
		rem.a = rem.a << (16 - binBytePos - sfc.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc.lb; binBytePos += sfc.lb; binBytePos &= 7;

		encode_block3(z->dat, z->len);
	}
	else if (wuns > (z->len) *f2)
	{
		SFC sfc = SFcode(3, 36);
		rem.a = sfc.code;
		rem.a = rem.a << (16 - binBytePos - sfc.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc.lb; binBytePos += sfc.lb; binBytePos &= 7;

		encode_block5(z->dat, z->len);

	}
	else
	{
		ESSC essc = encode_stationary_source_cbook(p);
		SFC sfc = SFcode(essc.codebook +4, 36);
		rem.a = sfc.code;
		rem.a = rem.a << (16 - binBytePos - sfc.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc.lb; binBytePos += sfc.lb; binBytePos &= 7;

		rem.a = lastBit;//编lastbit
		rem.a = rem.a << (15 - binBytePos);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr++; binBytePos++;
		binBytePos &= 7;

		encode_stationary_source_bin(runs->dat,runs->len, essc.codebook, 0, 0, 0, 0, 0);
	}
}