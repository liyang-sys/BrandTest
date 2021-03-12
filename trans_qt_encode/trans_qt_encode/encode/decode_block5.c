#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DEC decode_block5 (uchar *bin, int lenzbit)
{
   union data {
      uint a;
      uchar b[4];
   } rem;
   extern int ptr;
   int i, j;
   uchar x1, x2;
   int a; ///存放每次读bit的值
   int lenzbit5; /// lenzbit5是可用5block编码的z数组长度
   DEC dec1;

   dec1.z = (uchar *) calloc (lenzbit, sizeof(uchar));

   lenzbit5 = lenzbit - (lenzbit % 5);
   j = 0;
   x1 = ptr & 7;
   x2 = 0;
   while (j < lenzbit5) {
      a = (bin[ptr >> 3] >> (7 - x1)) & 1;
      if (a == 0) {
         ptr++; x1++; x1 &= 7;
         a = (bin[ptr >> 3] >> (7 - x1)) & 1;
         if (a == 0) {
            ptr++; x1++; x1 &= 7;
            a = (bin[ptr >> 3] >> (7 - x1)) & 1;
            if (a == 0) {
               /// 000: z = 00000;
               x2 += 5; x2 &= 7; j += 5;
            }
            else {
               ptr++; x1++; x1 &= 7;
               a = (bin[ptr >> 3] >> (7 - x1)) & 1;
               if (a == 0) {
                  /// 0010: z = 00001;
                  rem.a = 1;
                  rem.a = rem.a << (11 - x2);
                  dec1.z[j >> 3] |= rem.b[1];
                  dec1.z[(j >> 3) + 1] |= rem.b[0];
                  x2 += 5; x2 &= 7; j += 5;
               }
               else {
                  /// 0011: z = 00010;
                  rem.a = 2;
                  rem.a = rem.a << (11 - x2);
                  dec1.z[j >> 3] |= rem.b[1];
                  dec1.z[(j >> 3) + 1] |= rem.b[0];
                  x2 += 5; x2 &= 7; j += 5;
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
                  /// 0100: z = 00100;
                  rem.a = 4;
                  rem.a = rem.a << (11 - x2);
                  dec1.z[j >> 3] |= rem.b[1];
                  dec1.z[(j >> 3) + 1] |= rem.b[0];
                  x2 += 5; x2 &= 7; j += 5;
               }
               else {
                  /// 0101: z = 01000;
                  rem.a = 8;
                  rem.a = rem.a << (11 - x2);
                  dec1.z[j >> 3] |= rem.b[1];
                  dec1.z[(j >> 3) + 1] |= rem.b[0];
                  x2 += 5; x2 &= 7; j += 5;
               }
            }
            else {
               ptr++; x1++; x1 &= 7;
               a = (bin[ptr >> 3] >> (7 - x1)) & 1;
               if (a == 0) {
                  /// 0110: z = 10000;
                  rem.a = 16;
                  rem.a = rem.a << (11 - x2);
                  dec1.z[j >> 3] |= rem.b[1];
                  dec1.z[(j >> 3) + 1] |= rem.b[0];
                  x2 += 5; x2 &= 7; j += 5;
               }
               else {
                  ptr++; x1++; x1 &= 7;
                  a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                  if (a == 0) {
                     /// 01110: z = 00011;
                     rem.a = 3;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
                  else {
                     /// 01111: z = 00110;
                     rem.a = 6;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
               }
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
               ptr++; x1++; x1 &= 7;
               a = (bin[ptr >> 3] >> (7 - x1)) & 1;
               if (a == 0) {
                  ptr++; x1++; x1 &= 7;
                  a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                  if (a == 0)  {
                     /// 10000: z = 00101;
                     rem.a = 5;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
                  else {
                     /// 10001: z = 01100;
                     rem.a = 12;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
               }
               else {
                  ptr++; x1++; x1 &= 7;
                  a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                  if (a == 0) {
                     ///10010: z = 01010;
                     rem.a = 10;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
                  else {
                     /// 10011: z = 01001;
                     rem.a = 9;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
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
                     /// 10100: z = 11000;
                     rem.a = 24;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
                  else {
                     /// 10101: z = 10100;
                     rem.a = 20;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
               }
               else {
                  ptr++; x1++; x1 &= 7;
                  a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                  if (a == 0) {
                     /// 10110: z = 10010;
                     rem.a = 18;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
                  else {
                     /// 10111: z = 10001;
                     rem.a = 17;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
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
                  ptr++; x1++; x1 &= 7;
                  a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                  if (a == 0) {
                     /// 11000: z = 11100;
                     rem.a = 28;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
                  else {
                     /// 11001: z = 11001;
                     rem.a = 25;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
               }
               else {
                  ptr++; x1++; x1 &= 7;
                  a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                  if (a == 0) {
                     /// 11010: z = 11010;
                     rem.a = 26;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
                  }
                  else {
                     /// 11011: z = 00111;
                     rem.a = 7;
                     rem.a = rem.a << (11 - x2);
                     dec1.z[j >> 3] |= rem.b[1];
                     dec1.z[(j >> 3) + 1] |= rem.b[0];
                     x2 += 5; x2 &= 7; j += 5;
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
                     ptr++; x1++; x1 &= 7;
                     a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                     if (a == 0) {
                        /// 111000: z = 01011;
                        rem.a = 11;
                        rem.a = rem.a << (11 - x2);
                        dec1.z[j >> 3] |= rem.b[1];
                        dec1.z[(j >> 3) + 1] |= rem.b[0];
                        x2 += 5; x2 &= 7; j += 5;
                     }
                     else {
                        /// 111001: z = 01101;
                        rem.a = 13;
                        rem.a = rem.a << (11 - x2);
                        dec1.z[j >> 3] |= rem.b[1];
                        dec1.z[(j >> 3) + 1] |= rem.b[0];
                        x2 += 5; x2 &= 7; j += 5;
                     }
                  }
                  else {
                     ptr++; x1++; x1 &= 7;
                     a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                     if (a == 0) {
                        /// 111010: z = 01110;
                        rem.a = 14;
                        rem.a = rem.a << (11 - x2);
                        dec1.z[j >> 3] |= rem.b[1];
                        dec1.z[(j >> 3) + 1] |= rem.b[0];
                        x2 += 5; x2 &= 7; j += 5;
                     }
                     else {
                        /// 111011: z = 10011;
                        rem.a = 19;
                        rem.a = rem.a << (11 - x2);
                        dec1.z[j >> 3] |= rem.b[1];
                        dec1.z[(j >> 3) + 1] |= rem.b[0];
                        x2 += 5; x2 &= 7; j += 5;
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
                        /// 111100: z = 10110;
                        rem.a = 22;
                        rem.a = rem.a << (11 - x2);
                        dec1.z[j >> 3] |= rem.b[1];
                        dec1.z[(j >> 3) + 1] |= rem.b[0];
                        x2 += 5; x2 &= 7; j += 5;
                     }
                     else {
                        ptr++; x1++; x1 &= 7;
                        a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                        if (a == 0) {
                           /// 1111010: z = 10101;
                           rem.a = 21;
                           rem.a = rem.a << (11 - x2);
                           dec1.z[j >> 3] |= rem.b[1];
                           dec1.z[(j >> 3) + 1] |= rem.b[0];
                           x2 += 5; x2 &= 7; j += 5;
                        }
                        else {
                           /// 1111011: z = 11110;
                           rem.a = 30;
                           rem.a = rem.a << (11 - x2);
                           dec1.z[j >> 3] |= rem.b[1];
                           dec1.z[(j >> 3) + 1] |= rem.b[0];
                           x2 += 5; x2 &= 7; j += 5;
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
                           /// 1111100: z = 11101;
                           rem.a = 29;
                           rem.a = rem.a << (11 - x2);
                           dec1.z[j >> 3] |= rem.b[1];
                           dec1.z[(j >> 3) + 1] |= rem.b[0];
                           x2 += 5; x2 &= 7; j += 5;
                        }
                        else {
                           /// 1111101: z = 10111;
                           rem.a = 23;
                           rem.a = rem.a << (11 - x2);
                           dec1.z[j >> 3] |= rem.b[1];
                           dec1.z[(j >> 3) + 1] |= rem.b[0];
                           x2 += 5; x2 &= 7; j += 5;
                        }
                     }
                     else {
                        ptr++; x1++; x1 &= 7;
                        a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                        if (a == 0) {
                           /// 1111110: z = 01111;
                           rem.a = 15;
                           rem.a = rem.a << (11 - x2);
                           dec1.z[j >> 3] |= rem.b[1];
                           dec1.z[(j >> 3) + 1] |= rem.b[0];
                           x2 += 5; x2 &= 7; j += 5;
                        }
                        else {
                           ptr++; x1++; x1 &= 7;
                           a = (bin[ptr >> 3] >> (7 - x1)) & 1;
                           if (a == 0) {
                              /// 11111110: z = 11011;
                              rem.a = 27;
                              rem.a = rem.a << (11 - x2);
                              dec1.z[j >> 3] |= rem.b[1];
                              dec1.z[(j >> 3) + 1] |= rem.b[0];
                              x2 += 5; x2 &= 7; j += 5;
                           }
                           else {
                              /// 11111111: z = 11111;
                              rem.a = 31;
                              rem.a = rem.a << (11 - x2);
                              dec1.z[j >> 3] |= rem.b[1];
                              dec1.z[(j >> 3) + 1] |= rem.b[0];
                              x2 += 5; x2 &= 7; j += 5;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      ptr++; x1++; x1 &= 7;
   }
   /// 对未能5block编码的剩余z数组进行处理
   for (j = lenzbit5; j < lenzbit; j++) {
      a = (bin[ptr >> 3] >> (7 - x1)) & 1;
      ptr++; x1++; x1 &= 7;
      dec1.z[j >> 3] |= a << (7 - x2);
      x2++; x2 &= 7;
   }
   dec1.lenzbit = j;

   return dec1;
}
