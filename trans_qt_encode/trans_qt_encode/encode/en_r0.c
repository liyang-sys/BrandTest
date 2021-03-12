#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

void en_r0 (uint *r, int lenr)
{
    KIR kir1;
    SFC sfc1;
    EBNL ebnl1;
    GOLINV golinv1;
    ESSC essc1;
    union data {
      uint a;
      uchar b[4];
    } rem;
    extern int ptr;
    extern uchar *bin;
    int i ,x, x1;
    int lenzbit;
    float pw;
    float p0 = 0.4219;
    float p1 = 0.375;
    float p2 = 0.3125;
    int lastRun;
    int inv = 0;
    int flg;
    int nrs;
    int code;
    int lb;

    lenzbit = sum(r, lenr);

    pw = (float)lenr / (float)lenzbit;
    lastRun = r[lenr - 1];

    if (pw > 0.58)
    {
        inv = 1;
        kir1 = KW_inversion_runs (r, lenr);
        pw = (float)kir1.lenri / (float)lenzbit;
        r = kir1.ri;
        lenr = kir1.lenri;
    }

    if ((pw >= p0) && (inv == 0)) ///p0 = 0.4219
    {
        //printf("en_r0 1:\n");///
        golinv1 = GolombInv (r, 1, lenr);
        sfc1 = SFcode (1, 37);

        x = ptr & 7;
        bin[ptr >> 3] |= inv << (7 - x);
        ptr++; x++; x &= 7;
        rem.a = sfc1.code;
        rem.a = rem.a << (16 - x - sfc1.lb);
        bin[ptr >> 3] |= rem.b[1];
        bin[(ptr >> 3) + 1] |= rem.b[0];
        ptr += sfc1.lb; x += sfc1.lb; x &= 7;
        bin_add(bin, golinv1.z, golinv1.lenzbit);
    }
    else if (pw > p1) ///p1 = 0.375
    {
        if (inv == 0)
        {
            //printf("en_r0 2:\n");///
            golinv1 = GolombInv (r, 1, lenr);
            ebnl1 = encode_block_nrs_lenres (golinv1.z, golinv1.lenzbit, 3);
            sfc1 = SFcode (2, 37); lb = sfc1.lb; code = sfc1.code;
            sfc1 = SFcode (ebnl1.nrs + 1, 3); lb += sfc1.lb; code = (code << sfc1.lb) | sfc1.code;

            x = ptr & 7;
            bin[ptr >> 3] |= inv << (7 - x);
            ptr++; x++; x &= 7;
            rem.a = code;
            rem.a = rem.a << (16 - x - lb);
            bin[ptr >> 3] |= rem.b[1];
            bin[(ptr >> 3) + 1] |= rem.b[0];
            ptr += lb; x += lb; x &= 7;
            encode_block3 (golinv1.z, golinv1.lenzbit);
        }
        else
        {
            //printf("en_r0 3:\n");///
            golinv1 = GolombInv (r, 0, lenr);
            ebnl1 = encode_block_nrs_lenres (golinv1.z, golinv1.lenzbit, 3);
            sfc1 = SFcode (2, 37); lb = sfc1.lb; code = sfc1.code;
            ebnl1.nrs = ebnl1.lenres - ebnl1.nrs;
            sfc1 = SFcode (ebnl1.nrs + 1, 3); lb += sfc1.lb; code = (code << sfc1.lb) | sfc1.code;

            x = ptr & 7;
            bin[ptr >> 3] |= inv << (7 - x);
            ptr++; x++; x &= 7;
            rem.a = code;
            rem.a = rem.a << (16 - x - lb);
            bin[ptr >> 3] |= rem.b[1];
            bin[(ptr >> 3) + 1] |= rem.b[0];
            ptr += lb; x += lb; x &= 7;
            encode_block3 (golinv1.z, golinv1.lenzbit);
        }
    }
    else if (pw > p2) ///p2 = 0.3125
    {
        if (inv == 0)
        {
            //printf("en_r0 4:\n");///
            golinv1 = GolombInv (r, 1, lenr);
            ebnl1 = encode_block_nrs_lenres (golinv1.z, golinv1.lenzbit, 5);
            if ((ebnl1.lenres > 3) && (ebnl1.nrs == 1))
            {
                flg = 1;
            }
            else
            {
                flg = 0;
            }
            sfc1 = SFcode (3, 37); lb = sfc1.lb; code = sfc1.code;
            sfc1 = SFcode (ebnl1.nrs + 1, 5); lb += sfc1.lb + 1; code = (((code << sfc1.lb) | sfc1.code) << 1) + flg;

            x = ptr & 7;
            bin[ptr >> 3] |= inv << (7 - x);
            ptr++; x++; x &= 7;
            rem.a = code;
            rem.a = rem.a << (16 - x - lb);
            bin[ptr >> 3] |= rem.b[1];
            bin[(ptr >> 3) + 1] |= rem.b[0];
            ptr += lb; x += lb; x &= 7;
            encode_block5 (golinv1.z, golinv1.lenzbit);
        }
        else
        {
            //printf("en_r0 5:\n");///
            golinv1 = GolombInv (r, 0, lenr);
            ebnl1 = encode_block_nrs_lenres (golinv1.z, golinv1.lenzbit, 5);
            sfc1 = SFcode (3, 37); lb = sfc1.lb; code = sfc1.code;
            ebnl1.nrs = ebnl1.lenres - ebnl1.nrs;
            sfc1 = SFcode (ebnl1.nrs + 1, 5); lb += sfc1.lb; code = (code << sfc1.lb) | sfc1.code;

            x = ptr & 7;
            bin[ptr >> 3] |= inv << (7 - x);
            ptr++; x++; x &= 7;
            rem.a = code;
            rem.a = rem.a << (16 - x - lb);
            bin[ptr >> 3] |= rem.b[1];
            bin[(ptr >> 3) + 1] |= rem.b[0];
            ptr += lb; x += lb; x &= 7;
            encode_block5 (golinv1.z, golinv1.lenzbit);
        }
    }
    else
    {
        if ((lenr == 1) && (inv == 1))
        {
            //printf("ÒÑ¾­ÐÞ¸´bug en_r0 6:\n");///
            //bin[0] = 254;
            //ptr = 6;
			sfc1 = SFcode(37, 37);

			x = ptr & 7;
			bin[ptr >> 3] |= inv << (7 - x);
			ptr++; x++; x &= 7;
			rem.a = sfc1.code;
			rem.a = rem.a << (16 - x - sfc1.lb);
			bin[ptr >> 3] |= rem.b[1];
			bin[(ptr >> 3) + 1] |= rem.b[0];
			ptr += sfc1.lb; x += sfc1.lb; x &= 7;
        }
        else
        {
            essc1 = encode_stationary_source_cbook(pw);
            sfc1 = SFcode(essc1.codebook + 4, 37);

            x = ptr & 7;
            bin[ptr >> 3] |= inv << (7 - x);
            ptr++; x++; x &= 7;
            rem.a = sfc1.code;
            rem.a = rem.a << (16 - x - sfc1.lb);
            bin[ptr >> 3] |= rem.b[1];
            bin[(ptr >> 3) + 1] |= rem.b[0];
            ptr += sfc1.lb; x += sfc1.lb; x &= 7;
            encode_stationary_source_r_invr (r, inv, lastRun, pw, lenr, essc1.k, essc1.m, essc1.m1, essc1.m2, essc1.cls);
        }
    }
}
