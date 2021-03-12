#ifndef _TRIM_H
#define _TRIM_H

#define uchar unsigned char

int trim_coef(uchar *sgn, int *runs, int lg, int len, int ptr, float *fabsA);
int trim_coef_1pix(uchar *sgn, int *runs, int lg, int len, int ptr, float *fabsA);
int trim_coef_1pixf(uchar *sgn, int *runs, int lgr, int len, int ptr, float *fabsA);
int trim_coef_pix1(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA);
int trim_coef_2pix(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA);
int trim_coef_2pixf(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA);
int trim_coef_3pix(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA);
int trim_coef_3pixf(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA);
int trim_coef_4pix(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA);
int trim_coef_5pix(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA);
int trim_coef_pix1f(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA);
#endif