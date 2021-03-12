#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DEC decode_block3 (uchar *bin, int lenzbit)
{
   union data {
      unsigned short int a;
      uchar b[4];
   } rem;
   extern int ptr;
   int i, j;
   uchar x1, x2;
   int a; ///存放每次读bit的值
   int lenzbit3; ///lenzbit3是可用3block解码的z数组bit长度
   DEC dec1;

   dec1.z = (uchar *) calloc (lenzbit, sizeof(uchar));

   lenzbit3 = lenzbit - (lenzbit % 3);
   j = 0;
   x1 = ptr & 7;
   x2 = 0;
   while (j < lenzbit3) {
      a = (bin[ptr >> 3] >> (7 - x1)) & 1;
      if (a == 0) {
         ptr++; x1++; x1 &= 7;
         a = (bin[ptr >> 3] >> (7 - x1)) & 1;
         if (a == 0) {
            /// 00: z = 000;
            x2 += 3; x2 &= 7; j += 3;
         }
         else {
            ptr++; x1++; x1 &= 7;
            a = (bin[ptr >> 3] >> (7 - x1)) & 1;
            if (a == 0) {
               /// 010: z = 001;
               rem.a=1;
               rem.a = rem.a << (13 - x2);
               dec1.z[j >> 3] |= rem.b[1];
               dec1.z[(j >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; j += 3;
            }
            else {
               /// 011: z = 010;
               rem.a=2;
               rem.a = rem.a << (13 - x2);
               dec1.z[j >> 3] |= rem.b[1];
               dec1.z[(j >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; j += 3;
            }
         }
      }
      else {
         ptr++; x1++; x1 &= 7;
         a = (bin[ptr >> 3] >> (7 - x1)) & 1;
         if (a == 0) {
            ptr++; x1++; x1 &= 7;
            a = (bin[ptr >> 3] >> (7 - x1)) & 1;
            if (a == 0) {
               /// 100: z = 100;
               rem.a=4;
               rem.a = rem.a << (13 - x2);
               dec1.z[j >> 3] |= rem.b[1];
               dec1.z[(j >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; j += 3;
            }
            else {
               ///101: z = 110;
               rem.a=6;
               rem.a = rem.a << (13 - x2);
               dec1.z[j >> 3] |= rem.b[1];
               dec1.z[(j >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; j += 3;
            }
         }
         else {
            ptr++; x1++; x1 &= 7;
            a = (bin[ptr >> 3] >> (7 - x1)) & 1;
            if(a == 0){
               /// 110: z = 011;
               rem.a=3;
               rem.a = rem.a << (13 - x2);
               dec1.z[j >> 3] |= rem.b[1];
               dec1.z[(j >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; j += 3;
            }
            else {
               ptr++; x1++; x1 &= 7;
               a = (bin[ptr >> 3] >> (7 - x1)) & 1;
               if (a == 0) {
                  /// 1110: z = 101;
                  rem.a=5;
                  rem.a = rem.a << (13 - x2);
                  dec1.z[j >> 3] |= rem.b[1];
                  dec1.z[(j >> 3) + 1] |= rem.b[0];
                  x2 += 3; x2 &= 7; j += 3;
               }
               else {
                  /// 1111: z = 111;
                  rem.a=7;
                  rem.a = rem.a << (13 - x2);
                  dec1.z[j >> 3] |= rem.b[1];
                  dec1.z[(j >> 3) + 1] |= rem.b[0];
                  x2 += 3; x2 &= 7; j +=  3;
               }
            }
         }
      }
      ptr++; x1++; x1 &= 7;
   }
   /// 对未能3block编码的剩余z数组进行处理
   for (j = lenzbit3; j < lenzbit; j++) {
      a = (bin[ptr >> 3] >> (7 - x1)) & 1;
      ptr++; x1++; x1 &= 7;
      dec1.z[j >> 3] |= a << (7 - x2);
      x2++; x2 &= 7;
   }
   dec1.lenzbit = j;

   return dec1;
}
