#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

int bin_add (uchar *bin1, uchar *bin2, int lenbin2bit)
{
    union data {
      unsigned int a;
      uchar b[4];
    } rem;
    extern int ptr;
    int i, x, p;
    int lenbin2;
    x = ptr & 7;
    p = ptr;

    i = 0;
    lenbin2 = len_bit_to_byte(lenbin2bit);
    while (i < lenbin2)
    {
        rem.b[3] = bin2[i];
        rem.b[2] = bin2[i + 1];
        rem.b[1] = bin2[i + 2];
        rem.b[0] = bin2[i + 3];
        rem.a = rem.a >> x;
        bin1[p >> 3] |= rem.b[3];
        bin1[(p >> 3) + 1] |= rem.b[2];
        bin1[(p >> 3) + 2] |= rem.b[1];
        bin1[(p >> 3) + 3] |= rem.b[0];
        i += 3;
        p += 24;
    }
    ptr += lenbin2bit;
    x = ptr & 7;
    bin1[ptr >> 3] = (bin1[ptr >> 3] >> (7 - x)) << (7 - x); ///对多出来的bin置0
    bin1[(ptr >> 3) + 1] = 0;
    bin1[(ptr >> 3) + 2] = 0;
    bin1[(ptr >> 3) + 3] = 0;
}
