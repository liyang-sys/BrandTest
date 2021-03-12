#define _CRT_SECURE_NO_WARNINGS
#include "readYUV.h"
#include <stdio.h>
#include <stdlib.h>


unsigned char ***readYUV(char* filesourse, int h, int w, int f1, int Nframe, int fmt)
{
	FILE  *fp;
	int hh = h >> 1;
	int hw = w >> 1;
	int nby = 1;
	fp = fopen(filesourse, "rb");
	if (fp == NULL){
		printf("打开文件失败！\n");
		return NULL;
	}
	unsigned char*** ImgData = (unsigned char ***)malloc(sizeof(unsigned char **)* Nframe);
	for (int i = 0; i < Nframe; i++)
	{
		ImgData[i] = (unsigned char**)malloc(sizeof(unsigned char*)* h);
		for (int j = 0; j < h; j++)
		{
			ImgData[i][j] = (unsigned char*)malloc(sizeof(unsigned char)* w);
		}
	}

	for (int i = f1; i < f1 + Nframe; i++){
		if (fmt == 444){
			fseek(fp, (i - 1)*w*h * 3 * nby, 0);
			for (int j = 0; j < h; j++)
			{
				fread(&ImgData[i - f1][j][0], sizeof(unsigned char), w, fp);
			}
		}
		else if (fmt == 422){
			fseek(fp, (i - 1)*w*h * 2 * nby, 0);
			for (int j = 0; j < h; j++)
			{
				fread(&ImgData[i - f1][j][0], sizeof(unsigned char), w, fp);
			}
		}
		else{ //fmt=420
			fseek(fp, (i - 1)*w*h * 1.5 * nby, 0);
			for (int j = 0; j < h; j++)
			{
				fread(&ImgData[i - f1][j][0], sizeof(unsigned char), w, fp);
			}
		}
	}

	return ImgData;
}