#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

GOLINV GolombInv (uint *r, int lb, int lenr)
{
   GOLINV golinv1;
   int i;
   int x = 0;
   int lenzbit;
   golinv1.lenzbit = 0;

   lenzbit = sum(r, lenr);
   golinv1.z = (uchar *) calloc ((lenzbit >> 3) +1, sizeof(uchar));

   for (i = 0; i < (lenr - 1); i++) {
      while(r[i] > 1) {
         r[i]--;
         golinv1.lenzbit++;
         x++; x &= 7;
      }
      golinv1.z[golinv1.lenzbit >> 3] |= 1 << (7 - x);
      golinv1.lenzbit++; x++; x &= 7;
   }
   while(r[i] > 1) {
      r[i]--;
      golinv1.lenzbit++; x++; x &= 7;
   }
   golinv1.z[golinv1.lenzbit >> 3] |= lb << (7 - x);
   golinv1.lenzbit++;

   return golinv1;
}

