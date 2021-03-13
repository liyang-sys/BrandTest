#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

int ptr;
uchar *bin;


void test()
{
	ptr = 0;
	bin = calloc(409600, sizeof(uchar));
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
	FILE *fp = fopen("SN.txt", "rb");
	int snLen = 32400;
	Uint8_Dat sn;
	sn.len = snLen;
	sn.dat = (unsigned char *)calloc(snLen, sizeof(unsigned char));
	fread(sn.dat, sizeof(unsigned char), (snLen / 8) + 1, fp);
	fclose(fp);


	int cf0Len = 32400;
	float* cf0 = (float*)calloc(cf0Len, sizeof(float));
	fp = fopen("CF0Float.txt", "rb");
	fread(cf0, sizeof(int), cf0Len, fp);
	fclose(fp);
	int maxcf0 = 685;
	en_sub3d_sub2(cf0, &sn, cf0Len, maxcf0);
	writeBinToFile(bin, ptr);

	int lenbits = ptr;
	ptr = 0;
	int qctr = 0;
	DE_S_SUB sub = de_sub3d_sub2(bin, cf0Len, qctr, lenbits);


	int temp = 1;
}