#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

void en_KsideSub_n1s (uchar *z, int lenzbit, int lenr)
{
    union data {
      uint a;
      uchar b[4];
    } rem;
    extern int ptr;
    extern uchar *bin;
    int x;
    int i;
    float wuns;
    int flag;
    int lb;
    int code;
    float f0, f1, f2;
    SFC sfc1;
    GOL gol1;
    EBNL ebnl1;
    ESSC essc1;

    wuns = lenr;
    f0 = 0.421875, f1 = 0.375, f2 = 0.3125;

    if (wuns >= (lenzbit * f0))
    {
        //printf("en_KsideSub_n1s 1:\n");//
        ptr += 2;
        bin_add (bin, z, lenzbit);
    }
    else if (wuns > (lenzbit * f1))
    {
        //printf("en_KsideSub_n1s 2:\n");///
        ptr++;
        x = ptr & 7;
        bin[ptr >> 3] |= 1 << (7 - x);
        ptr++; x++; x &= 7;
        ebnl1 = encode_block_nrs_lenres (z, lenzbit, 3);
        sfc1 = SFcode(ebnl1.nrs + 1, 3); lb = sfc1.lb; code = sfc1.code;
        rem.a = code;
        rem.a = rem.a << (16 - x - lb);
        bin[ptr >> 3] |= rem.b[1];
        bin[(ptr >> 3) + 1] |= rem.b[0];
        ptr += lb; x += lb; x &= 7;
        encode_block3 (z, lenzbit);
    }
    else if (wuns > (lenzbit * f2))
    {
        //printf("en_KsideSub_n1s 3:\n");///
        x = ptr & 7;
        bin[ptr >> 3] |= 1 << (7 - x);
        ptr += 2; x += 2; x &= 7;
        ebnl1 = encode_block_nrs_lenres (z, lenzbit, 5);
        sfc1 = SFcode (ebnl1.nrs + 1, 5); lb = sfc1.lb; code = sfc1.code;
        rem.a = code;
        rem.a = rem.a << (16 - x - lb);
        bin[ptr >> 3] |= rem.b[1];
        bin[(ptr >> 3) + 1] |= rem.b[0];
        ptr += lb; x += lb; x &= 7;
        if ((ebnl1.lenres > 3) && (ebnl1.nrs == 1))
        {
            flag = 1;
        }
        else
        {
            flag = 0;
        }
        bin[ptr >> 3] |= flag << (7 - x);
        ptr++; x++; x &= 7;
        encode_block5 (z, lenzbit);
    }
    else if (wuns > (lenzbit * 0.242141716744801))
    {
        //printf("en_KsideSub_n1s 4:\n");///
        x = ptr & 7;
        bin[ptr >> 3] |= 1 << (7 - x);
        ptr++; x++; x &= 7;
        bin[ptr >> 3] |= 1 << (7 - x);
        ptr++; x++; x &= 7;
        gol1= Golomb0 (z, lenzbit);
        bin[ptr >> 3] |= gol1.lastBit << (7 - x);
        ptr++; x++; x &= 7;
        essc1 = encode_stationary_source_cbook (0.29);
        encode_stationary_source_bin (gol1.r, gol1.lenr, 0.29, essc1.k, essc1.m, essc1.m1, essc1.m2, essc1.cls);
    }
    else
    {
        //printf("en_KsideSub_n1s 5:\n");///
        //printf("p outside the range!\n");
    }
}

