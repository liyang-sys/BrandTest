#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"


void en_2n_code_sub(Uint32_Dat *cw, Uint32_Dat* ck,int n)
{
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
	extern int ptr;
	extern uchar *bin;
	int binBytePos = ptr & 7;

	if (n > 0.5)
	{
		if (cw != NULL)
		{
			for (int i = 0; i < cw->len; i++)
			{
				rem.a = cw->dat[i] - 1;
				rem.a = rem.a << (16 - n - binBytePos);
				bin[ptr >> 3] |= rem.b[1];
				bin[(ptr >> 3) + 1] |= rem.b[0];
				ptr += n; binBytePos += n;
				binBytePos &= 7;
			}
		}
		if (ck != NULL)
		{
			for (int i = 0; i < ck->len; i++)
			{
				rem.a = ck->dat[i] - 1;
				rem.a = rem.a << (16 - n - binBytePos);
				bin[ptr >> 3] |= rem.b[1];
				bin[(ptr >> 3) + 1] |= rem.b[0];
				ptr += n; binBytePos += n;
				binBytePos &= 7;
			}
		}
	}

}