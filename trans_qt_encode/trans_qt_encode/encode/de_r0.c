#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DER de_r0 (uchar *bin, int lenr, int lenbinbit)
{
    DER der1;
    DES des1;
    GOL gol1;
    DEC dec1;
    KIR kir1;
    extern int ptr;
    int i;
    int x, x1;
    int inv;
    int cbook;
    int lenc;
    int lencbit = 0;
    int n1s = 0;

    x = ptr & 7;
    inv = (bin[ptr >> 3] >> (7 - x)) & 1;
    ptr++; x++; x &= 7;
    des1 = deSFcode (bin, 37);
    cbook = des1.sym;

    if(cbook == 1)
    {
        //printf("decode 1:\n");//
        uchar *c = (uchar *) calloc ((100 * lenr), sizeof(uchar));
        x = ptr & 7;
        x1 = 0;
        while (n1s < lenr)
        {
            if (((bin[ptr >> 3] >> (7 - x)) & 1) > 0 )
            {
                c[lencbit >> 3] |= 1 << (7 - x1);
                n1s++;
            }
            ptr++; x++; x &= 7;
            lencbit++; x1++; x1 &= 7;
        }
        gol1 = Golomb0 (c, lencbit);
        der1.r = gol1.r;
        der1.lenr = gol1.lenr;
    }
    else if (cbook == 2)
    {
        if (inv == 0)
        {
            //printf("decode 2:\n");//
            dec1 = decode_block3_n1s (bin, lenr, lenbinbit);
            gol1 = Golomb0 (dec1.z, dec1.lenzbit);
        }
        else
        {
            //printf("decode 3:\n");//
            dec1 = decode_block3_n0s (bin, lenr, lenbinbit);
            gol1 = Golomb0 (dec1.z, dec1.lenzbit);
        }
        der1.r = gol1.r;
        der1.lenr = gol1.lenr;
    }
    else if (cbook == 3)
    {
        if (inv == 0)
        {
            //printf("decode 4:\n");//
            dec1 = decode_block5_n1s (bin, lenr, lenbinbit);
            gol1 = Golomb0 (dec1.z, dec1.lenzbit);
        }
        else
        {
            //printf("decode 5:\n");//
            dec1 = decode_block5_n0s (bin, lenr, lenbinbit);
            gol1 = Golomb0 (dec1.z, dec1.lenzbit);
        }
        der1.r = gol1.r;
        der1.lenr = gol1.lenr;
    }
    else
    {
        if ((cbook == 37) && (inv==1))
        {
            //printf("decode 6:\n");//
            der1.r = &lenr;
            der1.lenr = 1;
        }
        else
        {
            dec1 = decode_stationary_source_r_invr (bin, cbook - 4, lenr, inv, lenbinbit);
            der1.r = dec1.r;
            der1.lenr = dec1.lenr;
        }
    }

    if (inv == 1)
    {
        kir1 = KW_inversion_runs (der1.r, der1.lenr);
        der1.r = kir1.ri;
        der1.lenr = kir1.lenri;
    }

    return der1;
}
