#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

int en_runs_sep_sub (int thd, uint *rk, uint *nk, uint *rw, int lenrk, int lenrw)
{
    union data {
      uint a;
      uchar b[4];
    } rem;
    extern int ptr;
    extern uchar *bin;
    int x;
    int sumrk;
    float pk;
    SFC sfc1;
    GOLINV golinv1;

    en_r0 (rw, lenrw);

    if (thd > 1)
    {
        if ( (thd < 10 & lenrk > 500) | (thd < 7 & lenrk > 300) | (thd < 4 & lenrk > 150))
        {
            //printf("en_runs_sep_sub 1:\n");///
            en_Kside(rk, lenrk, thd);
        }
        else
        {
            sumrk = sum(rk, lenrk);
            pk = (float) lenrk / (float) sumrk;
            if (pk < 0.32)
            {
                //printf("en_runs_sep_sub 2:\n");///
                ptr++;
                encode_stationary_source_Nsym (rk, lenrk, pk, thd, 1, 16);
            }
            else
            {
                //printf("en_runs_sep_sub 3:\n");///
                x = ptr & 7;
                bin[ptr >> 3] |= 1 << (7 - x);
                ptr++;
                golinv1 = GolombInv (rk, 1, lenrk);
                en_KsideSub_n1s (golinv1.z, golinv1.lenzbit, lenrk);
            }
        }
    }
    else
    {
        //printf("en_runs_sep_sub 4:\n");///
    }

}
