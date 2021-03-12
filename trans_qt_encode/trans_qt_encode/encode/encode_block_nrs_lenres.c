#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

EBNL encode_block_nrs_lenres (uchar *z, int lenzbit, int j)
{
   EBNL ebnl1;
   int bit; ///取bit上的值
   int i, x;
   int lenzbit35; /// lenzbit35是可用3block或5block编码的z数组bit长度

   lenzbit35 = lenzbit - (lenzbit % j);
   ebnl1.lenres = lenzbit - lenzbit35;
   ebnl1.nrs = 0;

   /// 计算nrs、lenres
   x = lenzbit35 & 7;
   for (i = lenzbit35; i < lenzbit; i++) {
      bit  = (z[i >> 3] >> (7 - x)) & 1;
      ebnl1.nrs += bit;
      x++; x &= 7;
   }

   return ebnl1;
}
