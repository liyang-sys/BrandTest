#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

int len_bit_to_byte (int len_bit)
{
    int len_byte;
    if ((len_bit) % 8 == 0)
    {
        len_byte = len_bit >> 3;
    }
    else
    {
        len_byte = (len_bit >> 3) + 1;
    }
    return len_byte;
}
