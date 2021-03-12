#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"

#define uint unsigned int
#define uchar unsigned char

DEC decode_block5_n1s(uchar *bin, int n1s0, int lenbinbit)
{
	DEC dec1;
	DES des1;
	union data {
		uint a;
		uchar b[4];
	} rem;
	extern int ptr;
	uchar x1, x2;
	int bit; //存放每次读bit的值
	int nrs;
	int n1s = 0;
	int ia = 1;
	int n1ss;
	int ptr0;
	int ptr_inc;
	int nres;
	int flg;
	int temp1;
	int temp2;
	dec1.z = (uchar *)calloc(((5 * n1s0) >> 3) + 1, sizeof(uchar));

	des1 = deSFcode(bin, 5);
	nrs = des1.sym - 1;
	dec1.lenzbit = 0;
	x1 = ptr & 7;
	x2 = 0;
	flg = (bin[ptr >> 3] >> (7 - x1)) & 1;
	ptr++; x1++; x1 &= 7;
	temp1 = n1s0 - nrs;
	temp2 = 3 + nrs;
	/* 调试代码段 */
#if 0
	printf("  sym = %d\n", des1.sym);
	printf("  nrs = %d\n", nrs);
	printf("  flg = %d\n", flg);
	printf("  bin = 0x%x%x\n", bin[ptr >> 3], bin[1 + (ptr >> 3)]);
	printf("  ptr = %d\n", ptr);
	printf("  n1s0-nrs = %d\n", temp1); while (1);
#endif
	while ((n1s <= temp1) && (n1s < n1s0) && ((ptr + temp2) <= lenbinbit)) {
		ptr0 = ptr; n1ss = n1s;
		if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
			ptr++; x1++; x1 &= 7;
			if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
				ptr++; x1++; x1 &= 7;
				if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
					// 000: z = 00000;
					x2 += 5; x2 &= 7; dec1.lenzbit += 5;
				}
				else {
					ptr++; x1++; x1 &= 7;
					if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
						// 0010: z = 00001;
						rem.a = 1;
						rem.a = rem.a << (11 - x2);
						dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
						dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
						x2 += 5; x2 &= 7; dec1.lenzbit += 5;
						n1s += 1;
					}
					else {
						// 0011: z = 00010;
						rem.a = 2;
						rem.a = rem.a << (11 - x2);
						dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
						dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
						x2 += 5; x2 &= 7; dec1.lenzbit += 5;
						n1s += 1;
					}
				}
			}
			else {
				ptr++; x1++; x1 &= 7;
				if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
					ptr++; x1++; x1 &= 7;
					if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
						// 0100: z = 00100;
						rem.a = 4;
						rem.a = rem.a << (11 - x2);
						dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
						dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
						x2 += 5; x2 &= 7; dec1.lenzbit += 5;
						n1s += 1;
					}
					else {
						// 0101: z = 01000;
						rem.a = 8;
						rem.a = rem.a << (11 - x2);
						dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
						dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
						x2 += 5; x2 &= 7; dec1.lenzbit += 5;
						n1s += 1;
					}
				}
				else {
					ptr++; x1++; x1 &= 7;
					if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
						// 0110: z = 10000;
						rem.a = 16;
						rem.a = rem.a << (11 - x2);
						dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
						dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
						x2 += 5; x2 &= 7; dec1.lenzbit += 5;
						n1s += 1;
					}
					else {
						ptr++; x1++; x1 &= 7;
						if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
							// 01110: z = 00011;
							rem.a = 3;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
						else {
							// 01111: z = 00110;
							rem.a = 6;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
					}
				}
			}
		}
		else {
			ptr++; x1++; x1 &= 7;
			if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
				ptr++; x1++; x1 &= 7;
				if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
					ptr++; x1++; x1 &= 7;
					if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
						ptr++; x1++; x1 &= 7;
						if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
							// 10000: z = 00101;
							rem.a = 5;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
						else {
							// 10001: z = 01100;
							rem.a = 12;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
					}
					else {
						ptr++; x1++; x1 &= 7;
						if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
							// 10010: z = 01010;
							rem.a = 10;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
						else {
							// 10011: z = 01001;
							rem.a = 9;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
					}
				}
				else {
					ptr++; x1++; x1 &= 7;
					if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
						ptr++; x1++; x1 &= 7;
						if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
							// 10100: z = 11000;
							rem.a = 24;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
						else {
							// 10101: z = 10100;
							rem.a = 20;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
					}
					else {
						ptr++; x1++; x1 &= 7;
						if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
							// 10110: z = 10010;
							rem.a = 18;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
						else {
							// 10111: z = 10001;
							rem.a = 17;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 2;
						}
					}
				}
			}
			else {
				ptr++; x1++; x1 &= 7;
				if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
					ptr++; x1++; x1 &= 7;
					if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
						ptr++; x1++; x1 &= 7;
						if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
							// 11000: z = 11100;
							rem.a = 28;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 3;
						}
						else {
							// 11001: z = 11001;
							rem.a = 25;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 3;
						}
					}
					else {
						ptr++; x1++; x1 &= 7;
						if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
							// 11010: z = 11010;
							rem.a = 26;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 3;
						}
						else {
							// 11011: z = 00111;
							rem.a = 7;
							rem.a = rem.a << (11 - x2);
							dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
							dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
							x2 += 5; x2 &= 7; dec1.lenzbit += 5;
							n1s += 3;
						}
					}
				}
				else {
					ptr++; x1++; x1 &= 7;
					if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
						ptr++; x1++; x1 &= 7;
						if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
							ptr++; x1++; x1 &= 7;
							if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
								// 111000: z = 01011;
								rem.a = 11;
								rem.a = rem.a << (11 - x2);
								dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
								dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
								x2 += 5; x2 &= 7; dec1.lenzbit += 5;
								n1s += 3;
							}
							else {
								// 111001: z = 01101;
								rem.a = 13;
								rem.a = rem.a << (11 - x2);
								dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
								dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
								x2 += 5; x2 &= 7; dec1.lenzbit += 5;
								n1s += 3;
							}
						}
						else {
							ptr++; x1++; x1 &= 7;
							if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
								// 111010: z = 01110;
								rem.a = 14;
								rem.a = rem.a << (11 - x2);
								dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
								dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
								x2 += 5; x2 &= 7; dec1.lenzbit += 5;
								n1s += 3;
							}
							else {
								// 111011: z = 10011;
								rem.a = 19;
								rem.a = rem.a << (11 - x2);
								dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
								dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
								x2 += 5; x2 &= 7; dec1.lenzbit += 5;
								n1s += 3;
							}
						}
					}
					else {
						ptr++; x1++; x1 &= 7;
						if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
							ptr++; x1++; x1 &= 7;
							if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
								// 111100: z = 10110;
								rem.a = 22;
								rem.a = rem.a << (11 - x2);
								dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
								dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
								x2 += 5; x2 &= 7; dec1.lenzbit += 5;
								n1s += 3;
							}
							else {
								ptr++; x1++; x1 &= 7;
								if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
									// 1111010: z = 10101;
									rem.a = 21;
									rem.a = rem.a << (11 - x2);
									dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
									dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
									x2 += 5; x2 &= 7; dec1.lenzbit += 5;
									n1s += 3;
								}
								else {
									//1111011: z = 11110;
									rem.a = 30;
									rem.a = rem.a << (11 - x2);
									dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
									dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
									x2 += 5; x2 &= 7; dec1.lenzbit += 5;
									n1s += 4;
								}
							}
						}
						else {
							ptr++; x1++; x1 &= 7;
							if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
								ptr++; x1++; x1 &= 7;
								if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
									//1111100: z = 11101;
									rem.a = 29;
									rem.a = rem.a << (11 - x2);
									dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
									dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
									x2 += 5; x2 &= 7; dec1.lenzbit += 5;
									n1s += 4;
								}
								else {
									//1111101: z = 10111;
									rem.a = 23;
									rem.a = rem.a << (11 - x2);
									dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
									dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
									x2 += 5; x2 &= 7; dec1.lenzbit += 5;
									n1s += 4;
								}
							}
							else {
								ptr++; x1++; x1 &= 7;
								if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
									//1111110: z = 01111;
									rem.a = 15;
									rem.a = rem.a << (11 - x2);
									dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
									dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
									x2 += 5; x2 &= 7; dec1.lenzbit += 5;
									n1s += 4;
								}
								else {
									ptr++; x1++; x1 &= 7;
									if (((bin[ptr >> 3] >> (7 - x1)) & 1) == 0) {
										//11111110: z = 11011;
										rem.a = 27;
										rem.a = rem.a << (11 - x2);
										dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
										dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
										x2 += 5; x2 &= 7; dec1.lenzbit += 5;
										n1s += 4;
									}
									else {
										//11111111: z = 11111;
										rem.a = 31;
										rem.a = rem.a << (11 - x2);
										dec1.z[dec1.lenzbit >> 3] |= rem.b[1];
										dec1.z[(dec1.lenzbit >> 3) + 1] |= rem.b[0];
										x2 += 5; x2 &= 7; dec1.lenzbit += 5;
										n1s += 5;
									}
								}
							}
						}
					}
				}
			}
		}
		ptr++; x1++; x1 &= 7;
		ptr_inc = ptr - ptr0;
		ia++;
	}

	/* 调试代码段 */
#if 0
	for (int i = 0; i < (dec1.lenzbit + 7) / 8; i++) {
		printf("  (%3d):0x%2X", i + 1, dec1.z[i]);
		if ((i + 1) % 8 == 0) printf("\n");
	}
	printf("\n");
	printf("dec1.lenzbit = %d n1s0 = %d\n", dec1.lenzbit, n1s0);
	while (1);
#endif

	if ((flg == 1) || (ptr_inc == 3))
	{
		ia = ia - 2;
		ptr = ptr - 3; x1 = ptr & 7;
	}
	else
	{
		ia--;
	}

	if (nrs == 0)
	{
		if (n1s == n1s0)
		{
			dec1.lenzbit = ia * 5;
		}
		else
		{
			ia--;
			ptr = ptr - ptr_inc; x1 = ptr & 7;
			dec1.lenzbit = ia * 5;
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

		dec1.lenzbit = ia * 5; x2 = dec1.lenzbit & 7;
		while (n1s < n1s0)
		{
			bit = (bin[ptr >> 3] >> (7 - x1)) & 1;
			dec1.z[dec1.lenzbit >> 3] &= ~(1 << (7 - x2));
			dec1.z[dec1.lenzbit >> 3] |= bit << (7 - x2);
			//printf("bin = %d dec1.lenzbit = %d\n", ~(1<< (7-x2)), dec1.lenzbit + 1);
			n1s += bit;
			ptr++; x1++; x1 &= 7;
			dec1.lenzbit++; x2++; x2 &= 7;
		}
	}
	/* 调试代码段 */
#if 0
	for (int i = 0; i < (dec1.lenzbit + 7) / 8; i++) {
		printf("  (%3d):0x%2X", i + 1, dec1.z[i]);
		if ((i + 1) % 8 == 0) printf("\n");
	}
	printf("\n");
	printf("dec1.lenzbit = %d n1s0 = %d\n", dec1.lenzbit, n1s0);
	while (1);
#endif
	//printf("dec1.lenzbit = %d\n", dec1.lenzbit);
	return dec1;
}

