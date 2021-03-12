#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DEC de_KsideSub_n1s (uchar *bin, int n1s0, int lenbinbit)
{
    extern int ptr;
    int bn;
    int nc = 0;
    int n1s = 0;
    int x, x1;
    int lastBit;
    int lenc;
    int lencbit;
    //DEK dek1;
    DEC dec1;
    GOLINV golinv1;

    x = ptr & 7;
    bn = (bin[ptr >> 3] >> (7 - x)) & 1;
    ptr++; x++; x &= 7;
    bn = bn << 1;
    bn += (bin[ptr >> 3] >> (7 - x)) & 1;
    ptr++; x++; x &= 7;

    switch (bn)
    {
        case 0 :
        {
            //printf("de_KsideSub_n1s :  1\n");
            dec1.z = (uchar *) calloc (100*n1s0, sizeof(uchar));
            x1 = 0;
            while (n1s < n1s0)
            {
                if (((bin[ptr >> 3] >> (7 - x)) & 1) > 0)
                {
                    dec1.z[nc >> 3] |= 1 << (7 - x1);
                    n1s++;
                }
                ptr++; x++; x &= 7;
                nc++; x1++; x1 &= 7;
            }
            dec1.lenzbit = nc;
            break;
        }
        case 1 :
        {
            //printf("de_KsideSub_n1s :  2\n");
            dec1 = decode_block3_n1s (bin, n1s0, lenbinbit);
            break;
        }
        case 2 :
        {
            //printf("de_KsideSub_n1s :  3\n");
            dec1 = decode_block5_n1s (bin, n1s0, lenbinbit);
            break;
        }
        case 3 :
        {
            //printf("de_KsideSub_n1s :  4\n");
            lastBit = (bin[ptr >> 3] >> (7 - x));
            ptr++; x++; x &= 7;
            dec1 = decode_stationary_source_lenr (bin, 0, n1s0);
            golinv1 = GolombInv (dec1.r, lastBit, dec1.lenr);
            dec1.z = golinv1.z;
            dec1.lenzbit = golinv1.lenzbit;
            break;
        }
    }

    return dec1;
}

