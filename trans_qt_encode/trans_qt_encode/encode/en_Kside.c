#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

//sumr可以用nk来算吗
//int
void en_Kside(int *rk, int l_rk, int thd)
{
   extern int ptr;
   extern uchar *bin;
   int i, t0, t1, t;
   GOL g;
   ESSC e_cbook;
   SEP Separate0;
   int sumr = 0;
   float pk;
   int lens;                //中间值
   uchar codebook = 0;
   uchar *rA, *rB;    //把rk和rw转换为char的中间值


   switch (thd){
   case 2:
      rB = (uchar *)calloc((l_rk >> 3) + 1, sizeof(uchar));

      for (i = 0; i < l_rk; i++){
         rk[i]--;
         t = i & 7;
         rB[i >> 3] |= rk[i] << (7 - t);
         sumr += rk[i];
      }

      pk = (float)(sumr) / l_rk;

      if (pk>0.2421){
         ptr++;   //把0移到bin
         en_KsideSub(rB, l_rk, sumr);
      }
      else{
         g = Golomb0(rB, l_rk);
         t1 = ptr & 7;
         bin[ptr >> 3] |= (1 << (7 - t1));  //把1移到bin
         ptr++;

         //e_cbook = encode_stationary_source_cbook(pk, 1);
         e_cbook = encode_stationary_source_cbook(pk); ///修改

         if (e_cbook.codebook > 4){
            //printf("en_Kside 1:\n");///
            e_cbook.codebook = 4;
            for (i = 1; i >= 0; i--){
               t0 = ((e_cbook.codebook - 1) >> i) & 1;  //codebook的值直接赋给bin
               t1 = ptr & 7;
               bin[ptr >> 3] |= (t0 << (7 - t1));
               ptr++;
            }

            t1 = ptr & 7;
            bin[ptr >> 3] |= (g.lastBit << (7 - t1));  //把lastBit移到bin
            ptr++;

            encode_stationary_source_bin(g.r, g.lenr, codebook, 0, 0, 0, 0, 0);
         }
         else{
            //printf("en_Kside 2:\n");///
            e_cbook.codebook--;
            for (i = 1; i >= 0; i--){
               t0 = (e_cbook.codebook >> i) & 1;  //codebook的值直接赋给bin
               t1 = ptr & 7;
               bin[ptr >> 3] |= (t0 << (7 - t1));
               ptr++;
            }

            t1 = ptr & 7;
            bin[ptr >> 3] |= (g.lastBit << (7 - t1));  //把lastBit移到bin
            ptr++;

            encode_stationary_source_bin(g.r, g.lenr, pk, e_cbook.k, e_cbook.m, e_cbook.m1, e_cbook.m2, e_cbook.cls);   //bin
         }
         free(g.r);
      }
      free(rB);
      break;
   case 3:
      //printf("en_Kside 3:\n");///
      Separate0 = separate0(rk, 1, l_rk,0);
      en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);


      rA = (uchar *)calloc((Separate0.lrw >> 3) + 1, sizeof(uchar));
      sumr = 0;
      for (i = 0; i < Separate0.lrw; i++){
         Separate0.rw[i]--;
         t = i & 7;
         rA[i >> 3] |= Separate0.rw[i] << (7 - t);
         sumr += Separate0.rw[i];
      }

      en_KsideSub(rA, Separate0.lrw, sumr);

      free(Separate0.rw);
      free(Separate0.sep);
      free(rA);
      break;
   case 4:
      //printf("en_Kside 4:\n");///
      Separate0 = separate0(rk, 1, l_rk, 0);
      //enKsideSub =
      en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);
      lens = Separate0.lrw;

      Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
      //enKsideSub =
      en_KsideSub(Separate0.sep, lens, Separate0.lrw);

      rA = (uchar *)calloc((Separate0.lrw >> 3) + 1, sizeof(uchar));
      sumr = 0;
      for (i = 0; i < Separate0.lrw; i++){
         Separate0.rw[i]--;
         t = i & 7;
         rA[i >> 3] |= Separate0.rw[i] << (7 - t);
         sumr += Separate0.rw[i];
      }

      //enKsideSub =
      en_KsideSub(rA, Separate0.lrw, sumr);

      free(Separate0.rw);
      free(Separate0.sep);
      free(rA);
      break;
   case 5:
      //printf("en_Kside 5:\n");///
      Separate0 = separate0(rk, 1, l_rk, 0);

      //enKsideSub =
      en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);
      lens = Separate0.lrw;

      Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
      //enKsideSub =
      en_KsideSub(Separate0.sep, lens, Separate0.lrw);

      lens = Separate0.lrw;
      Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
      //enKsideSub =
      en_KsideSub(Separate0.sep, lens, Separate0.lrw);

      rA = (uchar *)calloc((Separate0.lrw >> 3) + 1, sizeof(uchar));
      sumr = 0;
      for (i = 0; i < Separate0.lrw; i++)
      {
         Separate0.rw[i]--;
         t = i & 7;
         rA[i >> 3] |= Separate0.rw[i] << (7 - t);
         sumr += Separate0.rw[i];
      }

      //enKsideSub =
      en_KsideSub(rA, Separate0.lrw, sumr);

      free(Separate0.rw);
      free(Separate0.sep);
      free(rA);
      break;
   case 6:
      //printf("en_Kside 6:\n");///
      Separate0 = separate0(rk, 2, l_rk,1);

      /*
      int s = 0;
      for (i = 0; i <= l_rk >> 3; i++)
         s += Separate0.sep[i];*/

      //enKsideSub =
      en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);
      lens = Separate0.lrw;

      rB = (uchar *)calloc((Separate0.lrk >> 3) + 1, sizeof(uchar));
      sumr = 0;
      for (i = 0; i < Separate0.lrk; i++){
         rk[i]--;
         t = i & 7;
         rB[i >> 3] |= rk[i] << (7 - t);
         sumr += rk[i];
      }

      //enKsideSub =
      en_KsideSub(rB, Separate0.lrk, sumr);

      Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
      //enKsideSub =
      en_KsideSub(Separate0.sep, lens, Separate0.lrw);
      lens = Separate0.lrw;

      Separate0 = separate0(Separate0.rw, 1, Separate0.lrw, 0);
      //enKsideSub =
      en_KsideSub(Separate0.sep, lens, Separate0.lrw);

      rA = (uchar *)calloc((Separate0.lrw >> 3) + 1, sizeof(uchar));
      sumr = 0;
      for (i = 0; i < Separate0.lrw; i++){
         Separate0.rw[i]--;
         t = i & 7;
         rA[i >> 3] |= Separate0.rw[i] << (7 - t);
         sumr += Separate0.rw[i];
      }

      //enKsideSub =
      en_KsideSub(rA, Separate0.lrw, sumr);

      free(Separate0.rw);
      free(Separate0.sep);
      free(rA);
      break;
   default:
      //printf("en_Kside 7:\n");///
      Separate0 = separate0(rk, 2, l_rk, 1);
      //enKsideSub =
      en_KsideSub(Separate0.sep, l_rk, Separate0.lrw);  //sumr要重新求

      rB = (uchar *)calloc((Separate0.lrk >> 3) + 1, sizeof(uchar));
      sumr = 0;
      for (i = 0; i < Separate0.lrk; i++){
         rk[i]--;
         t = i & 7;
         rB[i >> 3] |= rk[i] << (7 - t);
         sumr += rk[i];
      }

      //enKsideSub =
      en_KsideSub(rB, Separate0.lrk, sumr);

      sumr = 0;
      for (i = 0; i < Separate0.lrw; i++)
         sumr += Separate0.rw[i];

      //d =
      //encode_stationary_source_Nsym(Separate0.rw, Separate0.lrw, (float)Separate0.lrw / sumr, thd - 2, 1, 0);
	  encode_stationary_source_Nsym(Separate0.rw, Separate0.lrw, (float)Separate0.lrw / sumr, thd - 2, 0, 0);

      free(Separate0.rw);
      free(Separate0.sep);
      free(rB);
      break;
   }
}
