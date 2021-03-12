#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

void encode_block3 (uchar *z, int lenzbit)
{
   union data {
      unsigned short int a;
      uchar b[4];
   } rem;
   extern int ptr;
   extern uchar *bin;
   int i, x1, x2;//x1单指具体z的数据段单个字节上移到哪一位了,x2单指bin的数据段的单个字节移到哪一位了
   int bit; ///取bit上的值
   int lenzbit3; /// lenzbit3是可用3block编码的z数组bit长度

   lenzbit3 = lenzbit - (lenzbit % 3);
   x1 = 0;
   x2 = ptr & 7;
   for (i = 0; i < lenzbit3; i++)
   {
      if (((z[i >> 3] >> (7 - x1)) & 1) == 0)
      {
         i++; x1++; x1 &= 7;
         if (((z[i >> 3] >> (7 - x1)) & 1) == 0)
         {
            i++; x1++; x1 &= 7;
            if (((z[i >> 3] >> (7 - x1)) & 1) == 0)
            {
               ///000：bin = 00;
               ptr += 2; x2 += 2;
            }
            else
            {
               ///001：bin = 010;
               rem.a=2;
               rem.a = rem.a << (13 - x2);
               bin[ptr >> 3] |= rem.b[1];
               bin[(ptr >> 3) + 1] |= rem.b[0];
               ptr += 3; x2 += 3;
            }
         }
         else
         {
            i++; x1++; x1 &= 7;
            if (((z[i >> 3] >> (7 - x1)) & 1) == 0)
            {
               ///010：bin = 011;
               rem.a=3;
               rem.a = rem.a << (13 - x2);
               bin[ptr >> 3] |= rem.b[1];
               bin[(ptr >> 3) + 1] |= rem.b[0];
               ptr += 3; x2 += 3;
            }
            else
            {
               ///011：bin = 110;
               rem.a=6;
               rem.a = rem.a << (13 - x2);
               bin[ptr >> 3] |= rem.b[1];
               bin[(ptr >> 3) + 1] |= rem.b[0];
               ptr += 3; x2 += 3;
            }
         }
      }
      else
      {
         i++; x1++; x1 &= 7;
         if (((z[i >> 3] >> (7 - x1)) & 1) == 0)
         {
            i++; x1++; x1 &= 7;
            if (((z[i >> 3] >> (7 - x1)) & 1) == 0)
            {
               ///100：bin = 100;
               rem.a=4;
               rem.a = rem.a << (13 - x2);
               bin[ptr >> 3] |= rem.b[1];
               bin[(ptr >> 3) + 1] |= rem.b[0];
               ptr += 3; x2 += 3;
            }
            else
            {
               ///101：bin = 1110;
               rem.a=14;
               rem.a = rem.a << (12 - x2);
               bin[ptr >> 3] |= rem.b[1];
               bin[(ptr >> 3) + 1] |= rem.b[0];
               ptr +=  4; x2 += 4;
            }
         }
         else
         {
            i++; x1++; x1 &= 7;
            if (((z[i >> 3] >> (7 - x1)) & 1) == 0)
            {
               ///110：bin = 101;
               rem.a=5;
               rem.a = rem.a << (13 - x2);
               bin[ptr >> 3] |= rem.b[1];
               bin[(ptr >> 3) + 1] |= rem.b[0];
               ptr += 3; x2 += 3;
            }
            else
            {
               ///111：bin = 1111;
               rem.a=15;
               rem.a = rem.a << (12 - x2);
               bin[ptr >> 3] |= rem.b[1];
               bin[(ptr >> 3) + 1] |= rem.b[0];
               ptr +=  4; x2 += 4;
            }
         }
      }
      x1++; x1 &= 7;
      x2 &= 7;
   }

   /// 对未能3block编码的剩余z数组进行处理
   for (i = lenzbit3; i < lenzbit; i++)
   {
      bit  = (z[i >> 3] >> (7 - x1)) & 1;
      bin[ptr >> 3] |= bit << (7 - x2);
      x1++; x1 &= 7;
      ptr++; x2++; x2 &= 7;
   }
}
