#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DRSS de_runs_sep_sub (uchar *biny, int thd, int lenrk, int lenrw, int lenbinbit)
{
    extern int ptr;
    int i, x;
    int flg;
    DER der1;
    DES des1;
    //DEK dek1;
    DEC dec1;
    DEC dssnl1;
    DRSS drss1;
    GOL gol1;

    der1 = de_r0 (biny, lenrw, lenbinbit);
    drss1.rw = der1.r;
    drss1.lenrw = der1.lenr;

    if (thd > 1)
    {
        if (((thd < 10) && (lenrk > 500)) || ((thd < 7)  &&  (lenrk > 300)) || ((thd < 4) && (lenrk>150)))
        {
            //printf("de_runs_sep_sub 1:\n");
            dec1 = de_Kside (biny, lenrk, thd, lenbinbit);
            drss1.rk = dec1.r;
            drss1.lenrk = dec1.lenr;
        }
        else
        {
            x = ptr & 7;
            flg = ((biny[ptr >> 3] >> (7 - x)) & 1);
            ptr++;
            if (flg < 0.5)
            {
                //printf("de_runs_sep_sub 2:\n");
                des1 = deSFcode (biny, 16);
                des1.sym--; ///des1.sym = codebok
                dssnl1  = decode_stationary_source_Nsym_lenr (biny, des1.sym, lenrk, thd, lenbinbit);
                drss1.rk = dssnl1.r;
                drss1.lenrk = dssnl1.lenr;
            }
            else
            {
                //printf("de_runs_sep_sub 3:\n");
                dec1 = de_KsideSub_n1s (biny, lenrk, lenbinbit);
                gol1 = Golomb0 (dec1.z, dec1.lenzbit);
                drss1.rk = gol1.r;
                drss1.lenrk = gol1.lenr;
            }
        }
    }
    else
    {
        //printf("de_runs_sep_sub 4:\n");
        uint *rk = (uint *) malloc (lenbinbit * sizeof(uint));
        for (i = 0; i < lenrk; i++)
        {
            rk[i] = 1;
        }
        drss1.rk = rk;
        drss1.lenrk = lenrk;
    }

    return drss1;
}

