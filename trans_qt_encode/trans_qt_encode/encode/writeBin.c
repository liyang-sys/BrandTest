#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include "encoding.h"
#include "parameter_setting.h"
#include "all.h"

int ptr;
uchar *bin;
void writeBin()
{
	union data
	{
		unsigned short int a;
		uchar b[4];
	} rem;
	FILE *fp;
	ptr = 0;
	int binLen = 131072;
	bin = calloc((binLen / 8) + 1, sizeof(uchar));
	unsigned char *bigZ = (unsigned char *)calloc(binLen, sizeof(unsigned char));
	fp = fopen("de_zone_sub2B.txt", "rb");
	fread(bigZ, sizeof(unsigned char), binLen, fp);
	//for (int i = 0; i < binLen; i++)
	//{
	//	printf("%d", bigZ[i]);
	//}
	//printf("\n");
	int index = 0;
	int indexZ = 0;
	for (int i = 0; i < binLen; i++)
	{
		if (bigZ[i] == 1)
		{
			rem.a = 1;
			rem.a = rem.a << (15 - index);
			bin[indexZ >> 3] |= rem.b[1];
			bin[(indexZ >> 3) + 1] |= rem.b[0];
		}
		index++;
		indexZ++;
		index &= 7;
	}
	fclose(fp);

	fp = fopen("de_zone_sub2B_c.txt", "wb");
	fwrite(bin, sizeof(unsigned char), (binLen / 8) + 1, fp);
	fclose(fp);
}
