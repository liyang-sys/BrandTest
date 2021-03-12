#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

void en_KsideSub (uchar *c, int len, int wuns)  //len是c长度的8倍
{
   extern int ptr;
   extern uchar *bin;
   float f0 = 0.4219, f1 = 0.3750, f2 = 0.3125;
   int i, t, t0, t1;
   int inv = 0;
   GOL g;
   ESSC e_cbook;

   if ((wuns << 1)>len){
      for (i = 0; i <= (len >> 3); i++)
         c[i] = ~c[i];
      wuns = len - wuns;
      inv = 1;
   }
   if (wuns >= len*f0) {
      //printf("en_KsideSub 1:\n");///
      ptr += 2;
      for (i = 0; i < len; i++){   ////把c移到bin
         t1 = i & 7;
         t0 = ptr & 7;
         t = (c[i >> 3] >> (7 - t1)) & 1;   //取bin中的一位
         bin[ptr >> 3] |= (t << (7 - t0));
         ptr++;
      }

      t = ptr & 7;             //inv
      bin[ptr >> 3] |= inv << (7 - t);
      ptr++;
   }
   else if (wuns > len*f1){
      //printf("en_KsideSub 2:\n");///
      ptr++;         //01
      t = ptr & 7;
      bin[ptr >> 3] |= 1 << (7 - t);
      ptr++;

      //eb3 =
      encode_block3(c, len);

      /*
      for (i = 1265; i < 1275; i++)
         printf("%d ", bin[i]);
      int s = 0;
      for (i = 1265; i <= ptr >> 3; i++)
         s += bin[i];
*/
      t = ptr & 7;
      bin[ptr >> 3] |= inv << (7 - t);
      ptr++;
   }
   else if (wuns > len*f2){
      //printf("en_KsideSub 3:\n");///
      t = ptr & 7;  //10;
      bin[ptr >> 3] |= 1 << (7 - t);
      ptr += 2;

      //eb5 =
      encode_block5(c, len);

      t = ptr & 7;
      bin[ptr >> 3] |= inv << (7 - t);
      ptr++;

   }
   else if (wuns > len*0.2421){
      //printf("en_KsideSub 4:\n");///
      g = Golomb0(c, len);

      if (g.lastBit)
      for (i = 3; i >= 0; i--){    //把b0=13移到bin
         t1 = ptr & 7;
         t = (13 >> i) & 1;
         bin[ptr >> 3] |= (t << (7 - t1));
         ptr++;
      }
      else
      for (i = 3; i >= 0; i--){     //把b0=12移到bin
         t1 = ptr & 7;
         t = (12 >> i) & 1;
         bin[ptr >> 3] |= (t << (7 - t1));
         ptr++;
      }

      encode_stationary_source_bin(g.r, g.lenr, 0, 0, 0, 0, 0, 0);
	  //encode_stationary_source_bin(g.r, g.lenr, 0, 0, 2, 0, 0, 0);

      t = ptr & 7;
      bin[ptr >> 3] |= inv << (7 - t);
      ptr++;

      free(g.r);
   }
   else if (wuns > len*0.1797){
      //printf("en_KsideSub 5:\n");///
      g = Golomb0(c, len);

      if (g.lastBit)
      for (i = 4; i >= 0; i--){//把b0=13移到bin
         t1 = ptr & 7;
         t = (29 >> i) & 1;
         bin[ptr >> 3] |= (t << (7 - t1));
         ptr++;
      }
      else
      for (i = 4; i >= 0; i--){ //把b0=12移到bin
         t1 = ptr & 7;
         t = (28 >> i) & 1;
         bin[ptr >> 3] |= (t << (7 - t1));
         ptr++;
      }

     encode_stationary_source_bin(g.r, g.lenr, 1, 0, 0, 0, 0, 0);

      t = ptr & 7;
      bin[ptr >> 3] |= inv << (7 - t);
      ptr++;

      free(g.r);
    }
   else if (wuns > len*0.1294){
      //printf("en_KsideSub 6:\n");///
      g = Golomb0(c, len);

      if (g.lastBit)
      for (i = 5; i >= 0; i--)  {  //把b0=13移到bin
         t1 = ptr & 7;
         t = (61 >> i) & 1;
         bin[ptr >> 3] |= (t << (7 - t1));
         ptr++;
      }
      else
      for (i = 5; i >= 0; i--){    //把b0=12移到bin
         t1 = ptr & 7;
         t = (60 >> i) & 1;
         bin[ptr >> 3] |= (t << (7 - t1));
         ptr++;
      }

      encode_stationary_source_bin(g.r, g.lenr, 2, 0, 0, 0, 0, 0);//输入参数还有问题？

      t = ptr & 7;
      bin[ptr >> 3] |= inv << (7 - t);
      ptr++;

      free(g.r);
   }
   else{
      //printf("en_KsideSub 7:\n");///
      g = Golomb0(c, len);

      if (g.lastBit)
      for (i = 5; i >= 0; i--) { //把b0=13移到bin
         t1 = ptr & 7;
         t = (63 >> i) & 1;
         bin[ptr >> 3] |= (t << (7 - t1));
         ptr++;
      }
      else
      for (i = 5; i >= 0; i--){  //把b0=12移到bin
         t1 = ptr & 7;
         t = (62 >> i) & 1;
         bin[ptr >> 3] |= (t << (7 - t1));
         ptr++;
      }

      encode_stationary_source_bin(g.r, g.lenr, 3, 0, 0, 0, 0, 0);//输入参数还有问题？

      t = ptr & 7;
      bin[ptr >> 3] |= inv << (7 - t);
      ptr++;

      free(g.r);
    }

}
