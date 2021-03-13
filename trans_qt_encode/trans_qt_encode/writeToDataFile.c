#define _CRT_SECURE_NO_WARNINGS
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
void writeToDataFile(float* cf0, Uint8_Dat* sn, int lencf0)
{
	printf("cf0地址 = %p", cf0);
	FILE*fp = fopen("SN.txt", "wb");
	fwrite(sn->dat, sizeof(unsigned char), ((sn->len) / 8) + 1, fp);
	fclose(fp);
	int*cf0Int32 = calloc(lencf0, sizeof(int));
	//printf("64位int字节大小 = %d", sizeof(int));
	for (int i = 0; i < lencf0; i++)
	{
		cf0Int32[i] = cf0[i];
	}

	fp = fopen("CF0.txt", "wb");
	fwrite(cf0Int32, sizeof(int), lencf0, fp);
	fclose(fp);

	fp = fopen("CF0Float.txt", "wb");
	fwrite(cf0, sizeof(float), lencf0, fp);
	fclose(fp);

}