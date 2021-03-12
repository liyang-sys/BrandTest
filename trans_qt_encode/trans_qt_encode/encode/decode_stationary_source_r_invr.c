#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DEC decode_stationary_source_r_invr (uchar *bin, int codebook, int lenr, int inv, int lenbinbit)
{
    DEC dec1;

    if (inv == 0)
    {
        //printf("decode 7:\n");//
        dec1 = decode_stationary_source_lenr (bin, codebook, lenr);
    }
    else
    {
        //printf("decode 8&9:\n");//
        dec1 = decode_stationary_source_leninvr (bin, codebook, lenr, lenbinbit);
    }

    return dec1;
}
