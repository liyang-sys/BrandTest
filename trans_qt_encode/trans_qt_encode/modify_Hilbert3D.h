#ifndef _MODIFY_HIBERT3D_H_
#define _MODIFY_HIBERT3D_H_
#include<assert.h>

int nsubtem, L;
int W1;
int Hn_next, Wn_next, iw, ih;                //modify_Hilbert3D_ScanByBlock的返回值
int total_block;
int Nsub,N0;
typedef struct                               //存放分解后的序列、阶次、分解个数
{
	int block[20];
	int order[20];
	int num;
}block0;

block0 find_basic_block0(int length);						    //  find_2s_power
int Sign(int value);										//判断符号
block0 find_2s_inc0(int H0, int dH, block0 inc0);
int **create_arr0(int width, int height);                    //创建二维数组
int **draw_hibert(int W, int H);
int **modifyHilbert3D_BlockDeformation0(int rank, int wi, int high);
int *row0(int n);
void swap(int *a, int *b);
int **zigzag3D(int Hn, int dWn, int rank);
int *modify_Hilbert3D(int W, int H, int rank);
int **modify_Hilbert3D_sub(int dim, int inc, int **m0, int rank, int *st0);
int **modify_Hilbert3D_sub3D(int wi, int high, int dim, int inc, int **m0, int N0, int Nsub, int *st0, int rank);





#endif
