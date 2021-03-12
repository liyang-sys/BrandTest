//#include <stdio.h>
//#include "all.h"
//
//#define uint unsigned int
//#define uchar unsigned char
//
//GOL Golomb0 (uchar *z, int lenzbit)
//{
//   GOL gol1;
//   int i = 0;
//   int x = 0;
//   int n = 1; //计数
//
//   gol1.r = (uint *) calloc (lenzbit, sizeof(uint));
//
//   gol1.lenr = 0;
//   while (i < (lenzbit - 1)) {
//      if (((z[i >> 3] >> (7 - x)) & 1) == 0) {
//         n++;
//      }
//      else {
//         gol1.r[gol1.lenr] = n;
//         gol1.lenr++;
//         n = 1;
//      }
//      i++; x++; x &= 7; //x等于8时，归零
//   }
//   gol1.r[gol1.lenr] = n;
//   gol1.lenr++;
//
//   gol1.lastBit = z[i >> 3] >> (7 - x) & 1;
//   gol1.p = (float) gol1.lenr / (float) lenzbit;
//
//   return gol1;
//}
//
#include <stdio.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

GOL Golomb0(uchar *z, int lenzbit)
{
	GOL gol1;
	int i = 0;
	int x = 0;

	gol1.r = (uint *)calloc(lenzbit, sizeof(uint));

	gol1.lenr = 0;

	while (i < lenzbit) {
		gol1.r[gol1.lenr]++;
		if (((z[i >> 3] >> (7 - x)) & 1) == 1) {
			//printf("i = %d\n", i+1);
			gol1.lenr++;
		}
		i++; x++; x &= 7; //x等于8时，归零
	}
	gol1.lenr++;

	//gol1.lastBit = z[lenzbit >> 3] >> (7 - ((lenzbit - 1) & 7)) & 1;
	gol1.lastBit = z[(lenzbit - 1)>> 3] >> (7 - ((lenzbit - 1) & 7)) & 1; //修改的
	if (gol1.lastBit > 0)  gol1.lenr--;
	//printf("lastBit = %d\n", gol1.lastBit);

	gol1.p = (float)gol1.lenr / (float)lenzbit;

#if 0
	printf("r: len = %d\n", gol1.lenr);
	for (int i = 0; i < gol1.lenr; i++) {
		printf("(%3d):%4d ", i + 1, gol1.r[i]);
		if ((i + 1) % 8 == 0) printf("\n");
	}
	printf("\n");
	while (1);
#endif
	//printf("Nr = %d sum(runs) = %d gol1.lenr = %d p = %f\n", gol1.lenr ,lenzbit , gol1.lenr ,gol1.p ) ;while(1);

	return gol1;
}

