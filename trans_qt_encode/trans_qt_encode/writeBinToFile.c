#define _CRT_SECURE_NO_WARNINGS
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"
void writeBinToFile(unsigned char *bin,int ptr)
{
	FILE* fp = fopen("encodeZThd4.txt", "wb");
	fwrite(bin, sizeof(unsigned char), (ptr / 8) + 1, fp);
	fclose(fp);

}