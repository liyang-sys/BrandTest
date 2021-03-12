#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

//flg:whether or not calculate rk
SEP separate0(int *r, int thd, int len, int flg)
{
   int i, t;
   SEP Separate;

   Separate.rw = (unsigned int *)malloc(len*sizeof(unsigned int));
   Separate.sep = (uchar *)calloc((len >> 3) + 1, sizeof(uchar));

   if (flg)
   {
      Separate.lrk = 0, Separate.lrw = 0;
      for (i = 0; i < len; i++){
         if (r[i]>thd){

            t = i & 7;
            Separate.sep[i >> 3] |= (1 << (7 - t));  //把1移到sep

            Separate.rw[Separate.lrw] = r[i] - thd;          //把后面的减thd合并到这里
            Separate.lrw++;
         }
         else{
            r[Separate.lrk] = r[i];
            Separate.lrk++;
         }
      }
   }
   else
   {
      Separate.lrk = 0, Separate.lrw = 0;
      for (i = 0; i < len; i++){
         if (r[i]>thd){
            t = i & 7;
            Separate.sep[i >> 3] |= (1 << (7 - t));  //把1移到sep

            Separate.rw[Separate.lrw] = r[i] - thd;          //把后面的减thd合并到这里
            Separate.lrw++;
         }
      }
   }
   Separate.lensepbit = len;

   return Separate;
}
