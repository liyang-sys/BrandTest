#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

KIR KW_inversion_runs (uint *r, int lenr)
{
   KIR kir1;
   int i;
   int t = 1;
   kir1.lenri = 0;
   int lenzbit;

   lenzbit = sum(r, lenr);
   kir1.ri = (uint*) calloc ((lenzbit - lenr + 1), sizeof(uint));

   for (i = 0; i < lenr; i++)
   {
      if (r[i] > 1)
      {
         kir1.ri[kir1.lenri] = t;
         t = 1;
         while (r[i] > 1)
         {
             kir1.lenri++;
             kir1.ri[kir1.lenri] = 1;
             r[i]--;
         }
      }
      t++;
   }
   t--;
   kir1.ri[kir1.lenri] = t;
   kir1.lenri++;

   return kir1;
}
