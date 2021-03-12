#ifndef _QUANT_SUB_H_
#define _QUANT_SUB_H_

union Fabs{
	float **qcf;
	unsigned int **temp;
};

//*整理后的总的函数*//
void quant_subDC_harr(int a, int *w, int *h, int **m, float delta);
void quant_subDC_noharr(int a, int *w, int *h, int **m, float delta);
void quant_sub7_harr(int a, int *w, int *h, int **m, float delta);
void quant_sub7_noharr(int a, int *w, int *h, int **m, float delta);
void quant_sub8_harr(int a, int b, int *w, int *h, int **m, float delta);
void quant_sub8_noharr(int a, int b, int *w, int *h, int **m, float delta);

//*GetM*//获取路径 
int GetM(int W, int H, int *w, int *h, int **m);
int handleSn(unsigned char *sgn, int idx, int lg, int qctr, int *maxcf, float *fabsA);
void addSn(float *qcf, unsigned char *sgn, int lg);

#endif