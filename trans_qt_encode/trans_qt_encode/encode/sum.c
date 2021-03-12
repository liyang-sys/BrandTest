#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

int sum (uint *r, int lenr)
{
    int lenzbit = 0;
    int i;
    for (i = 0; i < lenr; i++)
    {
        lenzbit += r[i];
    }
    return lenzbit;
}
