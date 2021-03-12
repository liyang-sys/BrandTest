#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

void encode_stationary_source_r_invr (uint *r, int inv, int lastRun, float pw, int lenr, int k, int m, int m1, int m2, int cls)
{
    extern int ptr;
    extern uchar *bin;
    int i, x;

    if (inv == 0)
    {
        //printf("encode 7:\n");///
        encode_stationary_source_bin (r, lenr, pw, k, m, m1, m2, cls);
    }
    else
    {
        if (lastRun == 1)
        {
            //printf("encode 8:\n");///
            encode_stationary_source_bin (r, lenr, pw, k, m, m1, m2, cls);
            ptr++;
        }
        else
        {
            //printf("encode 9:\n");
            lenr = lenr - lastRun + 1;
            encode_stationary_source_bin (r, lenr, pw, k, m, m1, m2, cls);
            x = ptr & 7;
            for (i = 0; i < (lastRun - 1); i++){
                bin[ptr >> 3] |= 1 << (7 - x);
                ptr++; x++; x &= 7;
            }
            ptr++;
        }
    }
}
