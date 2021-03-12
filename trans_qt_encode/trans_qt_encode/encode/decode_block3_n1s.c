#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DEC decode_block3_n1s (uchar *bin, int n1s0, int lenbinbit)
{
   DEC dec1;
   DES des1;
   union data {
      unsigned short int a;
      uchar b[4];
   } rem;
   extern int ptr;
   int x1, x2;
   int bit; //存放每次读bit的值
   int nrs;
   int n1s = 0;
   int ia = 1;
   int n1ss;
   int ptr0;
   int ptr_inc;
   int nres;
   int temp1;
   int temp2;

   dec1.z = (uchar *) calloc (((5*n1s0) >> 3) + 1 , sizeof(uchar));

   des1 = deSFcode (bin, 3);
   nrs = des1.sym - 1;
   dec1.lenzbit = 0;
   x1 = ptr & 7;
   x2 = 0;
   temp1 = n1s0 - nrs;
   temp2 = nrs + 2;
   while ((n1s <= temp1) && (n1s < n1s0) && ((ptr + temp2) <= lenbinbit)) {
      ptr0 = ptr; n1ss = n1s;
      if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
         ptr++; x1++; x1 &= 7;
         if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
            // 00: z = 000;
            x2 += 3; x2 &= 7; dec1.lenzbit += 3;
         }
         else {
            ptr++; x1++; x1 &= 7;
            if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
               // 010: z = 001;
               rem.a=1;
               rem.a = rem.a << (13 - x2);
               dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
               dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; dec1.lenzbit += 3;
               n1s += 1;
            }
            else {
               // 011: z = 010;
               rem.a=2;
               rem.a = rem.a << (13 - x2);
               dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
               dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; dec1.lenzbit += 3;
               n1s += 1;
            }
         }
      }
      else {
         ptr++; x1++; x1 &= 7;
         if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
            ptr++; x1++; x1 &= 7;
            if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
               // 100: z = 100;
               rem.a=4;
               rem.a = rem.a << (13 - x2);
               dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
               dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; dec1.lenzbit += 3;
               n1s += 1;
            }
            else {
               //101: z = 110;
               rem.a=6;
               rem.a = rem.a << (13 - x2);
               dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
               dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; dec1.lenzbit += 3;
               n1s += 2;
            }
         }
         else {
            ptr++; x1++; x1 &= 7;
            if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
               // 110: z = 011;
               rem.a=3;
               rem.a = rem.a << (13 - x2);
               dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
               dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
               x2 += 3; x2 &= 7; dec1.lenzbit += 3;
               n1s += 2;
            }
            else {
               ptr++; x1++; x1 &= 7;
               if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
                  // 1110: z = 101;
                  rem.a=5;
                  rem.a = rem.a << (13 - x2);
                  dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
                  dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
                  x2 += 3; x2 &= 7; dec1.lenzbit += 3;
                  n1s += 2;
               }
               else {
                  // 1111: z = 111;
                  rem.a=7;
                  rem.a = rem.a << (13 - x2);
                  dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
                  dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
                  x2 += 3; x2 &= 7; dec1.lenzbit +=  3;
                  n1s += 3;
               }
            }
         }
      }
      ptr++; x1++; x1 &= 7;
      ptr_inc = ptr - ptr0;
      ia++;
   }
   ia--;

   if (nrs == 0)
   {
       if (n1s == n1s0)
       {
           dec1.lenzbit = ia * 3;
       }
       else
       {
           ia--;
           ptr = ptr - ptr_inc; x1 = ptr & 7;
           dec1.lenzbit = ia * 3;
       }
   }
   else
   {
       if (n1s > (n1s0 - nrs))
       {
           ia--;
           ptr = ptr - ptr_inc; x1 = ptr & 7;
           n1s = n1ss;
       }
       dec1.lenzbit = 3 * ia; x2 = dec1.lenzbit & 7;
       while (n1s < n1s0)
       {
           bit = (bin[ptr >> 3] >> (7 - x1)) & 1;
           dec1.z[dec1.lenzbit >> 3] |= bit << (7 - x2);
           n1s += bit;
           ptr++; x1++; x1 &= 7;
           dec1.lenzbit++; x2++; x2 &= 7;
       }
   }

   return dec1;
}

