#include "trim.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


static int len1;

int trim_coef(uchar *sgn, int *runs, int lg, int len, int ptr, float *fabsA)
{
	int lgr, lgl = 0, flg;
	lgr = trim_coef_1pix(sgn, runs, lg, len, ptr, fabsA);	//fabsA为quanTHD后的量化值
	lgl = trim_coef_1pixf(sgn, runs, lgr, len, ptr, fabsA);
	lgr = trim_coef_pix1(sgn, runs, lgl, lgr, len, ptr, fabsA);
	lgr = trim_coef_2pix(sgn, runs, lgl, lgr, len, ptr, fabsA);
	lgl = trim_coef_2pixf(sgn, runs, lgl, lgr, len, ptr, fabsA);
	lgr = trim_coef_3pix(sgn, runs, lgl, lgr, len, ptr, fabsA);
	lgl = trim_coef_3pixf(sgn, runs, lgl, lgr, len, ptr, fabsA);
	lgr = trim_coef_4pix(sgn, runs, lgl, lgr, len, ptr, fabsA);
	lgr = trim_coef_5pix(sgn, runs, lgl, lgr, len, ptr, fabsA);
	flg = trim_coef_pix1f(sgn, runs, lgl, lgr, len, ptr, fabsA);

	return flg;
}

int trim_coef_1pix(uchar *sgn, int *runs, int lg, int len, int ptr, float *fabsA)  //runs为原始的runs,
{
	int pos = ptr, i, lenp, j = 0;
	extern float *T;
	extern float *absA;
	int len1 = len + ptr - 1;
	pos += runs[0];
	pos--; //这里pos--是因为在matlab里，pos初始值是0；在matlab里pos+7，在这里相当于pos+6，因为数组索引是0开始
	for (i = 0; i < (lg - 1); i++){
		if (runs[i]>1){
			lenp = runs[i] < runs[i+1] ? runs[i] : runs[i+1]; //求出runs[i]和runs[i+1]最小值
			if (lenp>1){
				if (absA[pos] < T[lenp-1]){ //因为T和T2都已经增加了一个长度，所以这里和matlab一样，索引从1开始
					fabsA[pos] = 0;
					sgn[pos] = 0;
					pos += runs[i + 1];
					runs[i+1] += runs[i];
				}
				else{
					runs[j] = runs[i];
					pos += runs[i+1];
					j++;
				}
			}
			else{
				runs[j] = runs[i];
				pos += runs[i + 1];
				j++;
			}
		}
		else{
			runs[j] = runs[i];
			pos += runs[i + 1];
			j++;
		}
	}
	runs[j] = runs[lg-1];
	if (fabsA[len1] > 0){
		lenp = runs[lg - 1];
		if (lenp > 1){
			if (absA[len1] < T[lenp-1]){
				fabsA[len1] = 0;
				sgn[len1] = 0;
				runs[j]++;
			}
			else{
				j++;
				runs[j] = 1;
			}
		}
		else{
			j++;
			runs[j] = 1;
		}
	}
	else  runs[j]++;
	runs[0]--;
	return j;
}
/*************************************/
/*cf0和absA都fliplr, lgr为runs右边界,runs为从右往左开始,返回值为左边界lgl*/
int trim_coef_1pixf(uchar *sgn, int *runs, int lgr, int len, int ptr, float *fabsA) 
{
	int pos = len+ptr, i, lenp, j, st;
	extern float *T;
	extern float *absA;
	j = lgr;
	pos -= runs[lgr];
	if (!runs[0])
		st = 1;
	else st = 0;
	for (i = lgr; i > st; i--){
		if (runs[i]>1){
			lenp = runs[i] < runs[i-1] ? runs[i] : runs[i-1]; //求出runs[i]和runs[i-1]最小值
			if (lenp>1){
				if (absA[pos] < T[lenp-1]){
					fabsA[pos] = 0;
					sgn[pos] = 0;
					pos -= runs[i - 1];
					runs[ i - 1] += + runs[i];
				}
				else{
					runs[j] = runs[i];
					pos -= runs[i - 1];
					j--;
				}
			}
			else{
				runs[j] = runs[i];
				pos -= runs[i - 1];
				j--;
			}
		}
		else{
			runs[j] = runs[i];
			pos -= runs[i - 1];
			j--;
		}
	}
	runs[j] = runs[st];
	if (fabsA[ptr] > 0){
		lenp = runs[st];
		if (lenp > 1){
			if (absA[ptr] < T[lenp-1]){
				fabsA[ptr] = 0;
				sgn[ptr] = 0;
				runs[j]++;
			}
			else{
				j--;
				runs[j] = 1;
			}
		}
		else{
			j--;
			runs[j] = 1;
		}
	}
	else
		runs[j]++;
	runs[lgr]--;
	return j;
}
/*********************************************/
/*cf0和absA为正方向,lgl为runs左边界, lgr为右边界,runs为从左往右开始,返回值为lgr*/
int trim_coef_pix1(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA)  
{
	int pos = ptr, i, lenp, j, st;
	extern float *T;
	extern float *absA;
	int len1 = len + ptr - 1;
	j = lgl;
	pos += runs[lgl];
	pos--;
	if (!runs[lgr])  //判断右边界是否需要左移一位
		st = lgr - 1;
	else
		st = lgr;
	for (i = lgl; i < st; i++){
		if (runs[i]>1){
			lenp = runs[i]<runs[i + 1] ? runs[i] : runs[i + 1]; //求出runs[i]和runs[i+1]最小值
			if (lenp>1){
				if (absA[pos] < T[lenp-1]){
					fabsA[pos] = 0;
					sgn[pos] = 0;
					pos += runs[i + 1];
					runs[i + 1] += + runs[i];
				}
				else{
					runs[j] = runs[i];
					pos += runs[i + 1];
					j++;
				}
			}
			else{
				runs[j] = runs[i];
				pos += runs[i + 1];
				j++;
			}
		}
		else{
			runs[j] = runs[i];
			pos += runs[i + 1];
			j++;
		}
	}
	runs[j] = runs[st];
	if (fabsA[len1] > 0){
		lenp = runs[st];
		if (lenp > 1){
			if (absA[len1] < T[lenp-1]){
				fabsA[len1] = 0;
				sgn[len1] = 0;
			}
		}
	}
	return j;
}
/*********************************************/
/*cf0和absA为正方向,runs从左到右,返回的是右边界lgr*/
int trim_coef_2pix(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA)
{
	int i, j, pos = ptr, posA;
	int lenp1, lenp2, lenp3;
	int len1 = len + ptr - 1;
	extern float *T2;
	extern float *absA2;
	float ssb, th;
	i = lgl; j = lgl;
	pos += runs[i];
	pos--;
	while (i <= lgr-2 ){
		lenp1 = runs[i]; i++;
		lenp2 = runs[i]; i++;
		lenp3 = runs[i];
		if ((lenp1 > lenp2) && (lenp3 > lenp2)){
			posA = pos + lenp2;
			ssb = absA2[pos] + absA2[posA];
			if (T2[lenp1] < T2[lenp3])
				th = T2[lenp1] + T2[lenp2];
			else th = T2[lenp2] + T2[lenp3];
			if (ssb < th){
				fabsA[pos] = 0; fabsA[posA] = 0;
				sgn[pos] = 0; sgn[posA] = 0;
				runs[i] = lenp3 + lenp1 + lenp2;
				pos = posA + lenp3;
			}
			else{
				runs[j] = lenp1; j++;
				pos += lenp2; i--; //这里的i--是为了匹配前面的第一行和第二行的两个i++，因为前面相当于前进了两个，这里后退一个。
			}
		}
		else{
			runs[j] = lenp1; j++;
			pos += lenp2; i--;
		}
	}
	while (i <= lgr){
		runs[j] = runs[i];
		j++; i++;
	}
	if (!fabsA[len1]){
		j--;
		runs[j]++;
	}
	else
		runs[j] = 1;
	runs[lgl]--;
	return j;
}
/**********************************************/
/*cf0和absA都fliplr, lgr为runs右边界,runs为从右往左开始,返回值为左边界lgl*/
int trim_coef_2pixf(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA)
{
	int i, j, pos = len+ptr, posA, st;
	int lenp1, lenp2, lenp3;
	extern float *T2;
	extern float *absA2;
	float ssb, th;
	i = lgr; j = lgr;
	if (!runs[lgl]) //判断左边界是否右移一位
		st = lgl + 1;
	else
		st = lgl;
	pos -= runs[i];

	while (i >= st+2){
		lenp1 = runs[i]; i--;
		lenp2 = runs[i]; i--;
		lenp3 = runs[i];
		if ((lenp1 > lenp2) && (lenp3 > lenp2)){
			posA = pos - lenp2;
			ssb = absA2[pos] + absA2[posA];
			if (T2[lenp1] < T2[lenp3])
				th = T2[lenp1] + T2[lenp2];
			else  th = T2[lenp2] + T2[lenp3];
			if (ssb < th){
				fabsA[pos] = 0; fabsA[posA] = 0;
				sgn[pos] = 0; sgn[posA] = 0;
				runs[i] = lenp3 + lenp1 + lenp2;
				pos = posA - lenp3;
			}
			else{
				runs[j] = lenp1; j--;
				pos -= lenp2; i++;
			}
		}
		else{
			runs[j] = lenp1; j--;
			pos -= lenp2; i++;
		}
	}
	while (i >= st){
		runs[j] = runs[i];
		j--; i--;
	}
	if (!fabsA[ptr]){
		j++;
		runs[j]++;
	}
	else
		runs[j] = 1;
	runs[lgr]--;
	return j;
}
/*********************************************/
/*cf0和absA为正方向,runs从左到右,返回的是右边界lgr*/
int trim_coef_3pix(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA)
{
	int i, j, pos = ptr, posA, posB, lenp, st;
	int lenp1, lenp2, lenp3, lenp4;
	int len1 = len + ptr - 1;
	extern float *T2;
	extern float *absA2;
	float ssb, th;
	i = lgl; j = lgl;
	if (!runs[lgr])
		st = lgr - 1;
	else
		st = lgr;
	pos += runs[i];
	pos--;
	while (i <= st-3){
		lenp1 = runs[i]; i++;
		lenp2 = runs[i]; i++;
		lenp3 = runs[i]; i++;
		lenp4 = runs[i];
		lenp = lenp2 > lenp3 ? lenp2 : lenp3;
		if ((lenp1 > lenp) && (lenp4 > lenp)){
			posA = pos + lenp2; posB = posA + lenp3;
			ssb = absA2[pos] + absA2[posA] + absA2[posB];
			if (T2[lenp1] < T2[lenp4])
				th = T2[lenp1] + T2[lenp2] + T2[lenp3];
			else  th = T2[lenp2] + T2[lenp3] + T2[lenp4];
			if (ssb < th){
				fabsA[pos] = 0; fabsA[posA] = 0; fabsA[posB] = 0;
				sgn[pos] = 0; sgn[posA] = 0; sgn[posB] = 0;
				runs[i] = lenp4 + lenp1 + lenp2 + lenp3;
				pos = posB + lenp4;
			}
			else{
				runs[j] = lenp1; j++;
				pos += lenp2; i -= 2;
			}
		}
		else{
			runs[j] = lenp1; j++;
			pos += lenp2; i -= 2;
		}
	}
	while (i <= st){
		runs[j] = runs[i];
		j++; i++;
	}
	if (!fabsA[len1]){
		j--;
		runs[j]++;
	}
	else
		runs[j] = 1;
	runs[lgl]--;
	return j;
}
/************************************/
/*cf0和absA都fliplr, lgr为runs右边界,runs为从右往左开始,返回值为左边界lgl*/
int trim_coef_3pixf(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA)
{
	int i, j, pos=len+ptr, posA, posB, lenp, st;
	int lenp1, lenp2, lenp3, lenp4;
	extern float *T2;
	extern float *absA2;
	float ssb, th;
	i = lgr; j = lgr;
	if (!runs[lgl])
		st = lgl + 1;
	else
		st = lgl;
	pos -= runs[i];
	while (i >= st+3){
		lenp1 = runs[i]; i--;
		lenp2 = runs[i]; i--;
		lenp3 = runs[i]; i--;
		lenp4 = runs[i];
		lenp = lenp2 > lenp3 ? lenp2 : lenp3;
		if ((lenp1 > lenp) && (lenp4 > lenp)){
			posA = pos - lenp2; posB = posA - lenp3;
			ssb = absA2[pos] + absA2[posA] + absA2[posB];
			if (T2[lenp1] < T2[lenp4])
				th = T2[lenp1] + T2[lenp2] + T2[lenp3];
			else  th = T2[lenp2] + T2[lenp3] + T2[lenp4];
			if (ssb < th){
				fabsA[pos] = 0; sgn[pos] = 0;
				fabsA[posA] = 0; sgn[posA] = 0;
				fabsA[posB] = 0; sgn[posB] = 0;
				runs[i] = lenp4 + lenp1 + lenp2 + lenp3;
				pos = posB - lenp4;
			}
			else{
				runs[j] = lenp1; j--;
				pos -= lenp2; i += 2;
			}
		}
		else{
			runs[j] = lenp1; j--;
			pos -= lenp2; i += 2;
		}
	}
	while (i >= st){
		runs[j] = runs[i];
		j--; i--;
	}
	if (!fabsA[ptr]){
		j++;
		runs[j]++;
	}
	else
		runs[j] = 1;
	runs[lgr]--;
	return j;
}
/*********************************************/
/*cf0和absA为正方向,runs从左到右,返回的是右边界lgr*/
int trim_coef_4pix(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA)
{
	int i, j, pos = ptr, posA, posB, posC, lenp, st;
	int lenp1, lenp2, lenp3, lenp4, lenp5;
	extern float *T2;
	extern float *absA2;
	float ssb, th, temp1, temp2;
	i = lgl; j = lgl;
	if (!runs[lgr])
		st = lgr - 1;
	else
		st = lgr;
	pos += runs[i];
	pos--;
	while (i <= st-4){
		lenp1 = runs[i]; i++;
		lenp2 = runs[i]; i++;
		lenp3 = runs[i]; i++;
		lenp4 = runs[i]; i++;
		lenp5 = runs[i];
		if (lenp2 > lenp3) lenp = lenp2;
		else lenp = lenp3;
		if (lenp<lenp4) lenp = lenp4;
		if ((lenp1 > lenp) && (lenp5 > lenp)){
			posA = pos + lenp2; posB = posA + lenp3; posC = posB + lenp4;
			ssb = absA2[pos] + absA2[posA] + absA2[posB] + absA2[posC];
			temp1 = T2[lenp1] + T2[lenp2] + T2[lenp3] + T2[lenp4];
			temp2 = T2[lenp2] + T2[lenp3] + T2[lenp4] + T2[lenp5];
			if (temp1 < temp2) th = temp1;
			else th = temp2;
			if (ssb < th){
				fabsA[pos] = 0; sgn[pos] = 0;
				fabsA[posA] = 0; sgn[posA] = 0;
				fabsA[posB] = 0; sgn[posB] = 0;
				fabsA[posC] = 0; sgn[posC] = 0;
				runs[i] = lenp5 + lenp4 + lenp1 + lenp2 + lenp3;
				pos = posC + lenp5;
			}
			else{
				runs[j] = lenp1; j++;
				pos += lenp2; i -= 3;
			}
		}
		else{
			runs[j] = lenp1; j++;
			pos += lenp2; i -= 3;
		}
	}
	while (i <= st){
		runs[j] = runs[i]; 
		j++; i++;
	}
	j--;
	return j;
}
/*********************************************/
/*cf0和absA为正方向,runs从左到右,返回的是右边界lgr*/
int trim_coef_5pix(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA)
{
	int i, pos = ptr, j, posA, posB, posC, posD, lenp;
	int lenp1, lenp2, lenp3, lenp4, lenp5, lenp6;
	extern float *T2;
	extern float *absA2;
	float ssb, th, temp1, temp2;
	i = lgl; j = lgl;
	pos += runs[i];
	pos--;
	while (i <= lgr-5){
		lenp1 = runs[i]; i++;
		lenp2 = runs[i]; i++;
		lenp3 = runs[i]; i++;
		lenp4 = runs[i]; i++;
		lenp5 = runs[i]; i++;
		lenp6 = runs[i];
		if (lenp2 > lenp3) lenp = lenp2;
		else lenp = lenp3;
		if (lenp<lenp4) lenp = lenp4;
		if (lenp<lenp5) lenp = lenp5;
		if ((lenp1 > lenp) && (lenp6 > lenp)){
			posA = pos + lenp2; posB = posA + lenp3;
			posC = posB + lenp4; posD = posC + lenp5;
			ssb = absA2[pos] + absA2[posA] + absA2[posB] + absA2[posC] + absA2[posD];
			temp1 = T2[lenp1] + T2[lenp2] + T2[lenp3] + T2[lenp4] + T2[lenp5];
			temp2 = T2[lenp2] + T2[lenp3] + T2[lenp4] + T2[lenp5] + T2[lenp6];
			if (temp1 < temp2) th = temp1;
			else th = temp2;
			if (ssb < th){
				fabsA[pos] = 0; sgn[pos] = 0;
				fabsA[posA] = 0; sgn[posA] = 0;
				fabsA[posB] = 0; sgn[posB] = 0;
				fabsA[posC] = 0; sgn[posC] = 0;
				fabsA[posD] = 0; sgn[posD] = 0;
				runs[i] = lenp6 + lenp1 + lenp2 + lenp3 + lenp4 + lenp5;
				pos = posD + lenp6;
			}
			else{
				runs[j] = lenp1; j++;
				pos += lenp2; i -= 4;
			}
		}
		else{
			runs[j] = lenp1; j++;
			pos += lenp2; i -= 4;
		}
	}
	while (i <= lgr){
		runs[j] = runs[i];
		j++; i++;
	}
	j--;
	return j;
}
/*********************************************/
/*cf0和absA为正方向,runs从左到右,无返回值*/
int trim_coef_pix1f(uchar *sgn, int *runs, int lgl, int lgr, int len, int ptr, float *fabsA) //finally pix
{
	int pos = ptr, i, lenp, temp, j;
	int flg = 1;
	extern float *T;
	extern float *absA;
	int len1 = len + ptr - 1;
	pos += runs[lgl];
	pos--; j = 0;
	for (i = lgl; i <= lgr; i++){
		if (runs[i]>1){
			if (runs[i] < runs[i + 1]) lenp = runs[i];
			else lenp = runs[i + 1];
			if (lenp>1){
				if (absA[pos] < T[lenp-1]){
					fabsA[pos] = 0;
					sgn[pos] = 0;
					temp = runs[i + 1];
					runs[i + 1] = runs[i + 1] + runs[i];
					pos += temp;
				}
				else{
					pos += runs[i + 1];
					j++;
				}
			}
			else{
				pos += runs[i + 1];
				j++;
			}
		}
		else{
			pos += runs[i + 1];
			j++;
		}
	}
	if (fabsA[len1]>0){
		lenp = runs[lgr];
		if (lenp > 1){
			if (absA[len1] < T[lenp-1]){
				fabsA[len1] = 0;
				sgn[len1] = 0;
			}
		}

	}
	if (j < 2 && fabsA[len1] == 0){
		flg = 0;	//这里flg对应matlab里的max(cf0tr)>0的判断，j表示最后trim过后，runs的长度，
	}				//如果最后j=1，且最后一个值fabsA[len1]=0，那就说明整个子带都为零，即在matlab里max(cf0tr)==0
	free(runs);
	return flg;
}
