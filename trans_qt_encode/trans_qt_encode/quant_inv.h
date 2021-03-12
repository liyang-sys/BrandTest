#ifndef _QUANT_INV_H
#define _QUANT_INV_H

#define uchar unsigned char
#define uint unsigned int

float* rst_sub(uchar* bin, int* qcf, float delta, int lg);
//float* rst_sub(uchar* bin, float* qcf, float delta, int lg);
float* rstTHDctr(int *qcf, char *sgn, float T, float delta, int* nc, int maxcf, int lg);
float* rstEVENctr(int* qcf, char *sgn, float delta, int* nc, int maxcf, float ctr1, int lg);
void quant_DCinv(char* bin, int a, int lg, int* m, float delta);
void scan_inv(int a, int b, int c, int lg, int *m, int ptv, int harr, int dc, float *qcf);
int find_qctr(uchar *bin, int *qnt);

#endif