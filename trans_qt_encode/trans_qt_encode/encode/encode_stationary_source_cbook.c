#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

ESSC encode_stationary_source_cbook (float pw)
{
   ESSC essc1;
   essc1.codebook = 0, essc1.k = 0, essc1.m = 0, essc1.cls = 0;

   if (pw < 1 && pw != 0) {
      if (pw > 0.242141716744801)
         essc1.codebook = 0, essc1.k = 1, essc1.m = 2, essc1.cls = 0; //essc1.m=2
      else if (pw > 0.179664643992362)
         essc1.codebook = 1, essc1.k = 0, essc1.m = 3, essc1.m1 = 1, essc1.m2 = 2, essc1.cls = 1; //essc1.m=3
      else if (pw > 0.129449436703876)
         essc1.codebook = 2, essc1.k = 2, essc1.m = 4, essc1.cls = 0; //essc1.m=4
      else if (pw > 0.094276335736093)
         essc1.codebook = 3, essc1.k = 1, essc1.m = 6, essc1.m1 = 2, essc1.m2 = 4, essc1.cls = 1; //essc1.m=6

      else if (pw > 0.066967008463193)
         essc1.codebook = 4, essc1.k = 3, essc1.m = 8, essc1.cls = 0; //essc1.m=8
      else if (pw > 0.048304846989380)
         essc1.codebook = 5, essc1.k = 2, essc1.m = 12, essc1.m1 = 4, essc1.m2 = 8, essc1.cls = 1; //essc1.m=12
      else if (pw > 0.034063671075154)
         essc1.codebook = 6, essc1.k = 4, essc1.m = 16, essc1.cls = 0; //essc1.m=16
      else if (pw > 0.024451357947427)
         essc1.codebook = 7, essc1.k = 3, essc1.m = 24, essc1.m1 = 8, essc1.m2 = 16, essc1.cls = 1; //essc1.m=24

      else if (pw > 0.017179401454749)
         essc1.codebook = 8, essc1.k = 5, essc1.m = 32, essc1.cls = 0; //essc1.m=32
      else if (pw > 0.012301340462298)
         essc1.codebook = 9, essc1.k = 4, essc1.m = 48, essc1.m1 = 16, essc1.m2 = 32, essc1.cls = 1; //essc1.m=48
      else if (pw > 0.008626912537338)
         essc1.codebook = 10, essc1.k = 6, essc1.m = 64, essc1.cls = 0; //essc1.m=64
      else if (pw > 0.006169702847764)
         essc1.codebook = 11, essc1.k = 5, essc1.m = 96, essc1.m1 = 32, essc1.m2 = 64, essc1.cls = 1; //essc1.m=96
      else if (pw > 0.004322799566716)
         essc1.codebook = 12, essc1.k = 7, essc1.m = 128, essc1.cls = 0; //essc1.m=128
      else if (pw > 0.003089624313080)
         essc1.codebook = 13, essc1.k = 6, essc1.m = 192, essc1.m1 = 64, essc1.m2 = 128, essc1.cls = 1; //essc1.m=192
      else if (pw > 0.002163740670202)
         essc1.codebook = 14, essc1.k = 8, essc1.m = 256, essc1.cls = 0; //essc1.m=256
      else if (pw > 0.001546007225711)
         essc1.codebook = 15, essc1.k = 7, essc1.m = 384, essc1.m1 = 128, essc1.m2 = 256, essc1.cls = 1; //essc1.m=384

      else if (pw > 0.001082456190803)
         essc1.codebook = 16, essc1.k = 9, essc1.m = 512, essc1.cls = 0; //essc1.m=512
      else if (pw > 7.733e-4)
         essc1.codebook = 17, essc1.k = 8, essc1.m = 768, essc1.m1 = 256, essc1.m2 = 512, essc1.cls = 1; //essc1.m=768
      else if (pw > 5.4137e-4)
         essc1.codebook = 18, essc1.k = 10, essc1.m = 1024, essc1.cls = 0; //essc1.m=1024
      else if (pw > 3.8673e-4)
         essc1.codebook = 19, essc1.k = 9, essc1.m = 1536, essc1.m1 = 512, essc1.m2 = 1024, essc1.cls = 1; //essc1.m=1536
      else if (pw > 2.7072e-4)
         essc1.codebook = 20, essc1.k = 11, essc1.m = 2048, essc1.cls = 0; //essc1.m=2048
      else if (pw > 1.9338e-4)
         essc1.codebook = 21, essc1.k = 10, essc1.m = 3072, essc1.m1 = 1024, essc1.m2 = 2048, essc1.cls = 1; //essc1.m=3072
      else if (pw > 1.3537e-4)
         essc1.codebook = 22, essc1.k = 12, essc1.m = 4096, essc1.cls = 0; //essc1.m=4096
      else if (pw > 9.6696e-5)
         essc1.codebook = 23, essc1.k = 11, essc1.m = 6144, essc1.m1 = 2048, essc1.m2 = 4096, essc1.cls = 1; //essc1.m=6144

      else if (pw > 6.7688e-5)
         essc1.codebook = 24, essc1.k = 13, essc1.m = 8192, essc1.cls = 0; //essc1.m=8192
      else if (pw > 4.8349e-5)
         essc1.codebook = 25, essc1.k = 12, essc1.m = 12288, essc1.m1 = 4096, essc1.m2 = 8192, essc1.cls = 1; //essc1.m=12288
      else if (pw > 3.3845e-5)
         essc1.codebook = 26, essc1.k = 14, essc1.m = 16384, essc1.cls = 0; //essc1.m=16384
      else if (pw > 2.4175e-5)
         essc1.codebook = 27, essc1.k = 13, essc1.m = 24576, essc1.m1 = 8192, essc1.m2 = 16384, essc1.cls = 1; //essc1.m=24576
      else if (pw > 1.6922e-5)
         essc1.codebook = 28, essc1.k = 15, essc1.m = 32768, essc1.cls = 0; //essc1.m=32768
      else if (pw > 1.2087e-5)
         essc1.codebook = 29, essc1.k = 14, essc1.m = 49152, essc1.m1 = 16384, essc1.m2 = 32768, essc1.cls = 1; //essc1.m=49152
      else if (pw > 8.4612e-6)
         essc1.codebook = 30, essc1.k = 16, essc1.m = 65536, essc1.cls = 0; //essc1.m=65536
      else if (pw > 6.0437e-6)
         essc1.codebook = 31, essc1.k = 15, essc1.m = 98304, essc1.m1 = 32768, essc1.m2 = 65536, essc1.cls = 1; //essc1.m=98304
      else
         essc1.codebook = 32, essc1.k = 17, essc1.m = 131072, essc1.cls = 0; //essc1.m=131072
   }
   return essc1;
}
