#ifndef _QUANT_REST_H
#define _QUANT_REST_H

struct QuantTHD{
	int *runs;
	int lg;
};
typedef struct QuantTHD QT;

#define uint unsigned int
#define uchar unsigned char

int quanEVEN7(int *max, int a, int lg, int ptr, int harr, int *m, float delta0);
int quanEVEN8(int *max, int a, int b, int lg, int ptr, int harr, int *m, float delta0);
int quanEVEN2(int len, int ptr, float *absA, float delta);
QT quanTHD(uchar *sgn, int len, int ptr, float *absA, float delta);
int quanTHD2(float *absA, float delta);
int quanEVEN_DC(int a, int b, int c, int lg, int ptr, int ptv, int harr, int *m, float delta);
#endif