#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DEC decode_stationary_source_leninvr (uchar *bin, int codebook, int leninvr0, int lenbinbit)
{
   DEC dec1;
   union data3 {
      uint a;
      uchar b[4];
   } rem;
   int extern ptr;
   int i, x;
   int m1, m2, m, k;
   int a; //存放每次读bit的值
   int mai;
   int leninvr = 0;
   int flg;

   dec1.r = (uint *) calloc (lenbinbit, sizeof(uint));

   dec1.lenr = 0;
   x = ptr & 7;
   //m=2类型
   if (codebook << 31 == 0) {
      k = (codebook >> 1) + 1; m = 1 << k;
      while (leninvr < (leninvr0 - 1)) {
         mai = 0;
         while ((a = (bin[ptr >> 3] >> (7 - x)) & 1) == 1) {
            mai += m;
            x++; x &= 7; ptr++;
         }
         x++; x &= 7; ptr++;
         rem.a = 0;
         rem.b[2] |= bin[ptr >> 3]; rem.b[1] |= bin[(ptr >> 3) + 1]; rem.b[0] |= bin[(ptr >> 3) + 2];
         rem.a = rem.a << (8 + x); rem.a = rem.a >> (32 - k);
         dec1.r[dec1.lenr] = mai + rem.a + 1;
         leninvr += dec1.r[dec1.lenr] - 1;
         dec1.lenr++;
         x += k; x &= 7; ptr += k;
      }
   }
   //m=3类型
   else {
      k = codebook >> 1; m1 = 1 << k; m2 = m1 <<1; m = m1 + m2;
      while (leninvr < (leninvr0 - 1)) {
         mai = 0;
         while (((bin[ptr >> 3] >> (7 - x)) & 1) == 1) {
            mai += m;
            x++; x &= 7; ptr++;
         }
         x++; x &= 7; ptr++;
         if (k == 0) {
            if (((bin[ptr >> 3] >> (7 - x)) & 1) == 0) {
               dec1.r[dec1.lenr] = mai + 1;
            }
            else {
               x++; x &= 7; ptr++;
               if (((bin[ptr >> 3] >> (7 - x)) & 1) == 0) {
                  dec1.r[dec1.lenr] = mai + 2;
               }
               else {
                  dec1.r[dec1.lenr] = mai + 3;
               }
            }
            leninvr += dec1.r[dec1.lenr] - 1;//
            dec1.lenr++;
            x++; x &= 7; ptr++;
        }
        else {
           if (((bin[ptr >> 3] >> (7 - x)) & 1) ==0 ) {
              x++; x &= 7; ptr++;
              rem.a = 0;
              rem.b[2] |= bin[ptr >> 3]; rem.b[1] |= bin[(ptr >> 3) + 1]; rem.b[0] |= bin[(ptr >> 3) + 2];
              rem.a = rem.a << (8 + x); rem.a = rem.a >> (32 - k);
              dec1.r[dec1.lenr] = mai + rem.a + 1;
            }
            else {
               x++; x &= 7; ptr++;
               if (((bin[ptr >> 3] >> (7 - x)) & 1) == 0) {
                  x++; x &= 7; ptr++;
                  rem.a = 0;
                  rem.b[2] |= bin[ptr >> 3]; rem.b[1] |= bin[(ptr >> 3) + 1]; rem.b[0] |= bin[(ptr >> 3) + 2];
                  rem.a = rem.a << (8 + x); rem.a = rem.a >> (32 - k);
                  dec1.r[dec1.lenr] = mai + rem.a + 1 + m1;
               }
               else {
                  x++; x &= 7; ptr++;
                  rem.a = 0;
                  rem.b[2] |= bin[ptr >> 3]; rem.b[1] |= bin[(ptr >> 3) + 1]; rem.b[0] |= bin[(ptr >> 3) + 2];
                  rem.a = rem.a << (8 + x); rem.a = rem.a >> (32 - k);
                  dec1.r[dec1.lenr] = mai + rem.a + 1 + m2;
               }
            }
            leninvr += dec1.r[dec1.lenr] - 1;
            dec1.lenr++;
            x += k; x &= 7; ptr += k;
         }
      }
   }
   dec1.lenr--;

   flg = (bin[ptr >> 3] >> (7 - x)) & 1;
   ptr++; x++; x &= 7;
   if (flg == 1)
   {
       while (((bin[ptr >> 3] >> (7 - x)) & 1) == 1)
       {
           dec1.lenr++;
           dec1.r[dec1.lenr] = 1;
           ptr++; x++; x &= 7;
       }
       ptr++; x++; x &= 7;
       dec1.lenr++;
       dec1.r[dec1.lenr] = 1;
   }
   dec1.lenr++;

   return dec1;
}

