#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

void encode_block5 (uchar *z, int lenzbit)
{
   union data {
      uint a;
      uchar b[4];
   } rem;
   extern int ptr;
   extern uchar *bin;
   int i, x1, x2;
   int bit; ///取bit上的值
   int lenzbit5; /// lenzbit5是可用5block编码的z数组长度

   x1 = 0;
   x2 = ptr & 7;
   lenzbit5 = lenzbit - (lenzbit % 5);
   for (i = 0; i < lenzbit5; i++) {
      if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
         i++; x1++; x1 &= 7;
         if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
            i++; x1++; x1 &= 7;
            if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
               i++; x1++; x1 &= 7;
               if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///00000：bin = 000;
                     ptr += 3; x2 += 3;
                  }
                  else {
                     ///00001：bin = 0010;
                     rem.a = 2;
                     rem.a = rem.a << (12 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 4; x2 += 4;
                  }
               }
               else {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///00010：bin = 0011;
                     rem.a = 3;
                     rem.a = rem.a << (12 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 4; x2 += 4;
                  }
                  else {
                     ///00011：bin = 01110;
                     rem.a = 14;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
               }
            }
            else {
               i++; x1++; x1 &= 7;
               if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///00100：bin = 0100;
                     rem.a = 4;
                     rem.a = rem.a << (12 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 4; x2 += 4;
                  }
                  else {
                     ///00101：bin = 10000;
                     rem.a = 16;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
               }
               else {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///00110：bin = 01111;
                     rem.a = 15;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
                  else {
                     ///00111：bin = 11011;
                     rem.a = 27;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
               }
            }
         }
         else {
            i++; x1++; x1 &= 7;
            if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
               i++; x1++; x1 &= 7;
               if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///01000：bin = 0101;
                     rem.a = 5;
                     rem.a = rem.a << (12 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 4; x2 += 4;
                  }
                  else {
                     ///01001：bin = 10011;
                     rem.a = 19;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
               }
               else {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///01010：bin = 10010;
                     rem.a = 18;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
                  else {
                     ///01011：bin = 111000;
                     rem.a = 56;
                     rem.a = rem.a << (10 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 6; x2 += 6;
                  }
               }
            }
            else {
               i++; x1++; x1 &= 7;
               if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///01100：bin = 10001;
                     rem.a = 17;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
                  else {
                     ///01101：bin = 111001;
                     rem.a = 57;
                     rem.a = rem.a << (10 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 6; x2 += 6;
                  }
               }
               else {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///01110：bin = 111010;
                     rem.a = 58;
                     rem.a = rem.a << (10 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 6; x2 += 6;
                  }
                  else {
                     ///01111：bin = 1111110;
                     rem.a = 126;
                     rem.a = rem.a << (9 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 7; x2 += 7;
                  }
               }
            }
         }
      }
      else {
         i++; x1++; x1 &= 7;
         if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
            i++; x1++; x1 &= 7;
            if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
               i++; x1++; x1 &= 7;
               if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///10000：bin = 0110;
                     rem.a = 6;
                     rem.a = rem.a << (12 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 4; x2 += 4;
                  }
                  else {
                     ///10001：bin = 10111;
                     rem.a = 23;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
               }
               else {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///10010：bin = 10110;
                     rem.a = 22;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
                  else {
                     ///10011：bin = 111011;
                     rem.a = 59;
                     rem.a = rem.a << (10 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 6; x2 += 6;
                  }
               }
            }
            else {
               i++; x1++; x1 &= 7;
               if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///10100：bin = 10101;
                     rem.a = 21;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
                  else {
                     ///10101：bin = 1111010;
                     rem.a = 122;
                     rem.a = rem.a << (9 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 7; x2 += 7;
                  }
               }
               else {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///10110：bin = 111100;
                     rem.a = 60;
                     rem.a = rem.a << (10 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 6; x2 += 6;
                  }
                  else {
                     ///10111：bin = 1111101;
                     rem.a = 125;
                     rem.a = rem.a << (9 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 7; x2 += 7;
                  }
               }
            }
         }
         else {
            i++; x1++; x1 &= 7;
            if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
               i++; x1++; x1 &= 7;
               if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///11000：bin = 10100;
                     rem.a = 20;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
                  else {
                     ///11001：bin = 11001;
                     rem.a = 25;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
               }
               else {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///11010：bin = 11010;
                     rem.a = 26;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
                  else {
                     ///11011：bin = 11111110;
                     rem.a = 254;
                     rem.a = rem.a << (8 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 8; x2 += 8;
                  }
               }
            }
            else {
               i++; x1++; x1 &= 7;
               if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///11100：bin = 11000;
                     rem.a = 24;
                     rem.a = rem.a << (11 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 5; x2 += 5;
                  }
                  else {
                     ///11101：bin = 1111100;
                     rem.a = 124;
                     rem.a = rem.a << (9 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 7; x2 += 7;
                  }
               }
               else {
                  i++; x1++; x1 &= 7;
                  if (((z[i >> 3] >> (7 - x1)) & 1) == 0) {
                     ///11110：bin = 1111011;
                     rem.a = 123;
                     rem.a = rem.a << (9 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 7; x2 += 7;
                  }
                  else {
                     ///11111：bin = 11111111;
                     rem.a = 255;
                     rem.a = rem.a << (8 - x2);
                     bin[ptr >> 3] |= rem.b[1];
                     bin[(ptr >> 3) + 1] |= rem.b[0];
                     ptr += 8; x2 += 8;
                  }
               }
            }
         }
      }
      x1++; x1 &= 7;
      x2 &= 7;
   }

   /// 对未能5block编码的剩余z数组进行处理
   for (i = lenzbit5; i < lenzbit; i++) {
      bit  = (z[i >> 3] >> (7 - x1)) & 1;
      bin[ptr >> 3] |= bit << (7 - x2);
      x1++; x1 &= 7;
      ptr++; x2++; x2 &= 7;
   }
}

