#define _CRT_SECURE_NO_WARNINGS
#include "quant_sub.h"
#include "modify_Hilbert3D.h"
#include "modify_Hilbert2D.h"
#include <stdio.h>
#include <stdlib.h>
#include "quant_rest.h"
#include "trim.h"
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "quant_inv.h"
#include <trim.h>

#define NOENCODE 0
#define ENCODE 1
#define DECODE 1
#define INVQUANT 1
#define ulong unsigned long long
#define uchar unsigned char

void addSn(float *qcf, uchar *sgn, int lg)
{
	for (int i = 0; i < lg; i++){
		if (sgn[i] == 2)
			qcf[i] *= (-1);
	}

	return;
}

int quant_sub0(uchar *sgn, int max, int *nc, int lenth, int idx, float delta, float *fabsA)
{
	QT qtTHD; SFC se;
	uint tp = 1;
	int i, j, sum3 = 0, sum4 = 0, lenT = 131072, NT, lastlen, flg = 0, qctr, sumflg = 0;
	extern float *T;
	extern float *absA2;
	extern uchar *bin;
	extern int len;
	extern int ptr;
	float non, r61;
	float maxf0 = 0.0, temp, min, max_xx, sum1 = 0.0, sum2 = 0.0;
	float ctr[2], x2[14], x1[15], x[29], xx[30];

	if (max >= 6){
		int qnt = 1;
		if ((max < 1200) || (delta < 0.5)){  //因为这里的delta在外面已经进行了1/delta操作，所以是<0.5，而不>5
			non = 1.0 - (float)(nc[0] + nc[1] + nc[2] + nc[3] + nc[4]) / lenth;// non=sum(nc(6:lennc))/sum(nc);
			r61 = (float)nc[5] / nc[0];
			if (max >= 30){
				for (j = 0; j < 29; j += 2){
					x1[(j >> 1)] = nc[j] + nc[j + 1];
					if (j < 27){
						x2[(j >> 1)] = nc[j + 1] + nc[j + 2];
						x[j] = x1[(j >> 1)];
						x[j + 1] = x2[(j >> 1)];
					}
					else x[j] = x1[(j >> 1)];   //x[28]=x1[14]
				}
				min = x[0]; max_xx = 0;
				for (j = 0; j < 29; j++){	//找出x[]中的最小值
					if (min > x[j])
						min = x[j];
				}
				if (min > 0){
					for (j = 0; j < 28; j += 2){
						xx[j] = x[j + 1] / x[j];
						xx[j + 1] = x[j + 2] / x[j + 1];
					}
					for (j = 15; j < 25; j++){
						if (max_xx < xx[j]) max_xx = xx[j];  //求xx最大值
					}
					if (max_xx>0.9){
						for (j = 0; j < 30; j++)
							xx[j] = j + 1;
					}
				}
				else{
					for (j = 0; j < 30; j++)
						xx[j] = j + 1;
				}
			}
			else{
				for (j = 0; j < 30; j++)
					xx[j] = j + 1;
			}
			for (i = 15; i < 20; i++)//sum1为sum(xx(16:20)),sum2为sum(xx(21:25))
				sum1 += xx[i];
			for (i = 20; i < 25; i++)
				sum2 += xx[i];
			if (((non < 0.6) && (r61 < 0.71) && (sum1 < sum2)) || ((non < 0.36) && (r61 < 0.6))){
				qctr = 1;
				if (lenth>lenT){
					NT = (int)lenth / lenT;
					for (i = 0; i < NT; i++){
						qtTHD = quanTHD(sgn, lenT, idx, fabsA, delta);
						if ((qtTHD.lg > 1) || ((qtTHD.lg == 1) && (fabsA[lenT + idx - 1] != 0))){
							flg = trim_coef(sgn, qtTHD.runs, qtTHD.lg, lenT, idx, fabsA);
							sumflg += flg;
						}
						idx += lenT;
					}
					lastlen = lenth - NT*lenT;
					qtTHD = quanTHD(sgn, lastlen, idx, fabsA, delta);
					if ((qtTHD.lg > 1) || ((qtTHD.lg == 1) && (fabsA[lastlen + idx - 1] != 0))){
						flg = trim_coef(sgn, qtTHD.runs, qtTHD.lg, lastlen, idx, fabsA);
						sumflg += flg;
					}
				}
				else{
					qtTHD = quanTHD(sgn, lenth, idx, fabsA, delta);  //cf0=sb;
					if ((qtTHD.lg > 1) || ((qtTHD.lg == 1) && (fabsA[lenth + idx - 1] != 0))){
						flg = trim_coef(sgn, qtTHD.runs, qtTHD.lg, lenth, idx, fabsA);
						sumflg = flg;
					}
				}
				if (sumflg){ //sumflg大于0，表示max(cf0r)>0
					bin[ptr >> 3] = (3 << 6);
					ptr += 2;
				}
				else{
					bin[ptr >> 3] = 0;
					ptr += 1;
				}
			}
			else{
				for (i = 0; i < 5; i++)  sum3 += nc[i];
				for (i = 5; i < 15; i++) sum4 += nc[i];
				temp = (double)sum3 / sum4;
				if ((temp > 0.7071) || (((float)nc[0] / (float)nc[1] > 1.4) && (nc[0] > 35))){
					lenth += idx;
					for (i = idx; i < lenth; i++)
					{
						fabsA[i] = (int)(fabsA[i] * delta + 0.48); //quanTHD量化
						if (fabsA[i] == 0)
							sgn[i] = 0;
					}
					qctr = 1;
					bin[ptr >> 3] = (5 << 5);
					ptr += 3;
				}
				else{
					ctr[0] = quanEVEN2(lenth, idx, fabsA, delta);
					if (ctr[0]>63) ctr[0] = 63;
					if (ctr[0] < 0)  ctr[0] = 0;
					se = SFcode((uint)ctr[0]+1, 64);
					if (se.lb <= 5){
						bin[ptr >> 3] = (uchar)((4 << 5) | (se.code<<(5-se.lb)));
						ptr += (3 + se.lb);
					}
					else{
						uint temp = se.code << (5 + 8 - se.lb);
						uchar* ch = (uchar*)&temp;
						int id = ptr >> 3;
						bin[id] = ((uchar)(4 << 5) | ch[1]);
						bin[id + 1] = ch[0];
						ptr += (3 + se.lb);
					}
					qctr = 0;
				}
			}
		}
		else{
			for (i = 0; i < 5; i++)  sum3 += nc[i];
			for (i = 5; i < 15; i++) sum4 += nc[i];
			temp = (double)sum3 / sum4;
			if (((temp > 0.7071) && max < 2600) || (((float)nc[0] / (float)nc[1] > 1.4) && nc[0]>35)){
				lenth += idx;
				for (i = idx; i < lenth; i++)
				{
					fabsA[i] = (int)(fabsA[i] * delta + 0.48); //quantTHD量化,相比较matlab,通过计算省去了offset
					if (!fabsA[i])
						sgn[i] = 0;
				}
				qctr = 1;
				bin[ptr >> 3] = (5 << 5);
				ptr += 3;
			}
			else{
				ctr[0] = quanEVEN2(lenth, idx, fabsA, delta);
				if (ctr[0]>63) ctr[0] = 63;
				if (ctr[0] < 0)  ctr[0] = 0;
				se = SFcode((uint)ctr[0]+1, 64);
				if (se.lb <= 5){
					bin[ptr >> 3] = (uchar)((4 << 5) | (se.code << (5 - se.lb)));
					ptr += (3 + se.lb);
				}
				else{
					uint temp = se.code << (5 + 8 - se.lb);
					uchar* ch = (uchar*)&temp;
					bin[ptr >> 3] = (4 << 5) | ch[1];
					bin[(ptr >> 3) + 1] = ch[0];
					ptr += (3 + se.lb);
				}
				qctr = 0;
			}
		}
	}
	else{
		qctr = 1;
		bin[ptr >> 3] = 0; ptr += 1;
	}
	memset(nc, 0, (int)(len)*sizeof(int));
	return qctr;
}

/*
编码DC系数，返回值表示编码总位数
*/
int quant_DC(int a, int lg, int *m, float delta)
{
	int b = 0, c = 0, i = 0, qf;
	ulong max, min, temp;
	extern uchar *bin;
	extern float ***PTVData;
	extern union Fabs f1;
	int offset0 = 8;
	int offset1 = 32;
	max = 0; min = PTVData[a][b][c];
	for (int j = 0; j < lg; j++){
		f1.qcf[0][j] = PTVData[a][b][c];
		f1.temp[0][j] &= 0x7fffffff; //将符号位置零求绝对值
		qf = (int)(f1.qcf[0][j] * delta + 0.5); //round的量化
		f1.qcf[0][j] = qf;
		if (qf > max)
			max = qf;
		if (qf < min)
			min = qf;
		switch (m[i])
		{
		case 3:
			c += offset0;
			break;
		case -3:
			c -= offset0;
			break;
		case 2:
			b += offset0;
			break;
		case -2:
			b -= offset0;
			break;
		case 1:
			a += offset1;
			break;
		case -1:
			a -= offset1;
			break;
		default:
			break;
		}i++;
	}
	int nbits, nbitsMax, Nsym;
	int lnbMin, lnbMax, resMin, resMax; //lnbMin,lnbMax表示占了多少字节，resMin,resMax表示最后一个字节占了多少位
	int idx = 0, reslb;
	nbits = (int)(17.5 + log(delta) / log(2.0)); //这里delta=1/delta，所以减号变加号
	nbitsMax = (int)(log(max) / log(2)) + 1;
	//编码Maxqsb
	resMax = nbits % 8;
	if (resMax){
		lnbMax = (int)(nbits / 8 + 1);
		temp = max << (8 - resMax);
		unsigned char *ch = (unsigned char*)&temp;
		for (int i = lnbMax - 1; i >= 0; i--){
			bin[idx] = ch[i]; 
			idx++;
		}
	}
	else{
		lnbMax = (int)(nbits / 8);
		temp = max;
		unsigned char *ch = (unsigned char*)&temp;
		for (int i = lnbMax - 1; i >= 0; i--){
			bin[idx] = ch[i];
			idx++;
		}
	}
	//编码Minqsb
	resMin = nbitsMax % 8;
	if (resMin){
		lnbMin = (int)(nbitsMax / 8 + 1);
		if (resMax){
			temp = min << (16 - resMax - resMin);
			unsigned char *ch = (unsigned char*)&temp;
			bin[idx - 1] |= ch[lnbMin];
			if (resMin <= (8 - resMax)){
				for (int i = lnbMin-1; i > 0; i--){
					bin[idx] = ch[i];
					idx++;
				}
				reslb = 8 - resMax - resMin;  //reslb表示编码后一个字节还空多少位
			}
			else{
				for (int i = lnbMin-1; i >= 0; i--){
					bin[idx] = ch[i];
					idx++;
				}
				reslb = 16 - resMin - resMax;
			}
		}
		else{
			temp = min << (8 - resMin);
			unsigned char *ch = (unsigned char*)&temp;
			for (int i = lnbMin - 1; i >= 0; i--){
				bin[idx] = ch[i];
				idx++;
			}
			reslb = 8 - resMin;
		}
	}
	else{
		lnbMin = (int)(nbitsMax / 8);
		if (resMax){
			temp = min << (8 - resMax);
			unsigned char *ch = (unsigned char*)&temp;
			bin[idx - 1] |= ch[lnbMin];
			for (int i = lnbMin-1; i >= 0; i--){
				bin[idx] = ch[i];
				idx++;
			}
			reslb = resMax;
		}
		else{
			temp = min;
			unsigned char *ch = (unsigned char*)&temp;
			for (int i = lnbMin - 1; i >= 0; i--){
				bin[idx] = ch[i];
				idx++;
			}
			reslb = 0;
		}
	}
	//编码系数
	int lenb, resbit;  //lenb表示占了多少个字节，resbit表示最后一个字节用了多少位
	SFC se;
	Nsym = max - min + 1;
	for (int ii = 0; ii < lg; ii++){
		qf = f1.qcf[0][ii] + 1 - min;
		se = SFcode(qf, Nsym);
		temp = se.code;
		resbit = se.lb % 8;
		if (reslb){  //1、上一次的编码最后一个字节还有reslb个空位
			if (resbit){ //1.1、此系数最后一个字节需要用到resbit个位
				temp <<= (reslb + 8 - resbit); //左移，方便将上一次编码的空位填满
				unsigned char *ch = (unsigned char*)&temp;
				lenb = (int)(se.lb >> 3) + 1;
				bin[idx - 1] |= ch[lenb];//通过“或操作”将上一次编码的空位填满
				if (resbit <= reslb){	//需要用到的位数小于或等于空位
					for (int j = lenb - 1; j > 0; j--){
						bin[idx] = ch[j];	//将剩余位依次填入下一个字节
						idx++;
					}
					reslb = reslb - resbit;	//计算本次编码后，最后一个字节还剩多少个空位
				}
				else{	//需要的位数大于空位
					for (int j = lenb - 1; j >= 0; j--){ //原理同上
						bin[idx] = ch[j];
						idx++;
					}
					reslb = 8 - (resbit - reslb);
				}
			}
			else{	//1.2、此系数刚好需要整数个字节
				temp <<= reslb; //原理同上
				unsigned char *ch = (unsigned char*)&temp;
				lenb = (int)(se.lb >> 3);
				bin[idx - 1] |= ch[lenb];
				for (int j = lenb - 1; j >= 0; j--){
					bin[idx] = ch[j];
					idx++;
				}
				reslb = reslb;
			}
		}
		else{	//2、最后一个字节没有空位
			if (resbit){	//2.1、此系数最后一个字节需要用到resbit个位
				temp <<= (8 - resbit);
				unsigned char *ch = (unsigned char*)&temp;
				lenb = (int)(se.lb >> 3) + 1;
				for (int j = lenb - 1; j >= 0; j--){
					bin[idx] = ch[j];
					idx++;
				}
				reslb = 8 - resbit;
			}
			else{	//2.2、此系数刚好需要整数个字节
				unsigned char *ch = (unsigned char*)&temp;
				lenb = (int)(se.lb >> 3);
				for (int j = lenb - 1; j >= 0; j--){
					bin[idx] = ch[j];
					idx++;
				}
				reslb = 0;
			}
		}
	}
	return idx * 8 - reslb;
}

void quant_subDC_harr(int a, int *w, int *h, int **m, float delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern float ***rePTVData;
	extern int **nc;
	extern int len;
	extern union Fabs f1;
	int idx = 0, max, lg, qctr, maxcf0 = 0;
	float delta0;
	Uint8_Dat sign;

	DE_S_SUB desub;
	int lenbits;
	int qnt = 0;
	float *qcf0;

	delta0 = (float)1.0 / delta;

	/*第三次harr变换后的低频部分*/
//DC--------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	lg = w[0] * h[0];//第三次PTV，对应matlab里Ldc=coef(:,:,2)
	int lenb = quant_DC(a, lg, m[0], delta0);
#if NOENCODE
	ptr = 0;
	quant_DCinv(bin, a, lg, m[0], delta);
#endif
#if INVQUANT
	ptr = 0;
	quant_DCinv(bin, a, lg, m[0], delta);
#endif

//2---------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar)); //注意将bin清零，否则影响后面结果
	idx = 0; ptr = 0; //ptr清零
	lg = w[1] * h[0];
	max = quanEVEN_DC(a, 0, 4, lg, idx, 3, 3, m[1], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]); //将符号sn用位来储存，0表示正，1表示负，返回值为总位数
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);
		free(desub.sn.dat); free(desub.cfw.dat);
		//FILE *fp2;
		//fp2 = fopen("E:\\程序代码\\量化\\2021-01-15\\test_en_sub3d_sub2_sub\\decf0.txt", "wb");
		//fwrite(desub.cf.dat, sizeof(int), desub.cf.len, fp2);/*4*/
		//fclose(fp2);
		//fp2 = fopen("E:\\程序代码\\量化\\2021-01-15\\test_en_sub3d_sub2_sub\\desn.txt", "wb");
		//fwrite(desub.sn.dat, sizeof(int), desub.sn.len, fp2);/*4*/
		//fclose(fp2);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, 0, 4, lg, m[1], 3, 3, 1, qcf0);	//反向扫描
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 4, lg, m[1], 3, 3, 1, qcf0);	//反向扫描
#endif
//3--------------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[1] * h[1];
	max = quanEVEN_DC(a, 4, 4, lg, idx, 3, 3, m[3], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);	//量化
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);	//编码
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a, 4, 4, lg, m[3], 3, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a, 4, 4, lg, m[3], 3, 3, 1, qcf0);	//反向扫描
#endif
//4-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[0] * h[1];
	max = quanEVEN_DC(a, 4, 0, lg, idx, 3, 3, m[2], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a, 4, 0, lg, m[2], 3, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a, 4, 0, lg, m[2], 3, 3, 1, qcf0);	//反向扫描
#endif
//5--------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[3] * h[2];//第二次PTV，对应matlab里L5=coef(:,:,5)
	max = quanEVEN_DC(a, 0, 2, lg, idx, 2, 3, m[4], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a, 0, 2, lg, m[4], 2, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a, 0, 2, lg, m[4], 2, 3, 1, qcf0);	//反向扫描
#endif	
//6-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[2] * h[3];
	max = quanEVEN_DC(a, 2, 0, lg, idx, 2, 3, m[5], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a, 2, 0, lg, m[5], 2, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a, 2, 0, lg, m[5], 2, 3, 1, qcf0);	//反向扫描
#endif
//7--------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[3] * h[3];
	max = quanEVEN_DC(a, 2, 2, lg, idx, 2, 3, m[6], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a, 2, 2, lg, m[6], 2, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a, 2, 2, lg, m[6], 2, 3, 1, qcf0);	//反向扫描
#endif
//8-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[5] * h[4];//第一次PTV
	max = quanEVEN_DC(a, 0, 1, lg, idx, 1, 3, m[7], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a, 0, 1, lg, m[7], 1, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a, 0, 1, lg, m[7], 1, 3, 1, qcf0);	//反向扫描
#endif
//9---------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[4] * h[5];
	max = quanEVEN_DC(a, 1, 0, lg, idx, 1, 3, m[8], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a, 1, 0, lg, m[8], 1, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a, 1, 0, lg, m[8], 1, 3, 1, qcf0);	//反向扫描
#endif
//10---------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[5] * h[5];
	max = quanEVEN_DC(a, 1, 1, lg, idx, 1, 3, m[9], delta0*10.0);	
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a, 1, 1, lg, m[9], 1, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a, 1, 1, lg, m[9], 1, 3, 1, qcf0);	//反向扫描
#endif

	/*第三次harr变换后的高频部分*/
//1-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[0] * h[0];//第三次PTV
	max = quanEVEN_DC(a + 16, 0, 0, lg, idx, 3, 3, m[0], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 16, 0, 0, lg, m[0], 3, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 0, 0, lg, m[0], 3, 3, 1, qcf0);	//反向扫描
#endif
//2-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[1] * h[0];
	max = quanEVEN_DC(a + 16, 0, 4, lg, idx, 3, 3, m[1], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 16, 0, 4, lg, m[1], 3, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 0, 4, lg, m[1], 3, 3, 1, qcf0);	//反向扫描
#endif
//3-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[1] * h[1];
	max = quanEVEN_DC(a + 16, 4, 4, lg, idx, 3, 3, m[3], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 16, 4, 4, lg, m[3], 3, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 4, 4, lg, m[3], 3, 3, 1, qcf0);	//反向扫描
#endif
//4---------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[0] * h[1];
	max = quanEVEN_DC(a + 16, 4, 0, lg, idx, 3, 3, m[2], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a+16, 4, 0, lg, m[2], 3, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 4, 0, lg, m[2], 3, 3, 1, qcf0);	//反向扫描
#endif
//5--------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[3] * h[2];//第二次PTV
	max = quanEVEN_DC(a + 16, 0, 2, lg, idx, 2, 3, m[4], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a+16, 0, 2, lg, m[4], 2, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 0, 2, lg, m[4], 2, 3, 1, qcf0);	//反向扫描
#endif
//6-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[2] * h[3];
	max = quanEVEN_DC(a + 16, 2, 0, lg, idx, 2, 3, m[5], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 16, 2, 0, lg, m[5], 2, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 2, 0, lg, m[5], 2, 3, 1, qcf0);	//反向扫描
#endif
//7-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0; 
	lg = w[3] * h[3];
	max = quanEVEN_DC(a + 16, 2, 2, lg, idx, 2, 3, m[6], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 16, 2, 2, lg, m[6], 2, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 2, 2, lg, m[6], 2, 3, 1, qcf0);	//反向扫描
#endif
//8-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[5] * h[4];//第一次PTV
	max = quanEVEN_DC(a + 16, 0, 1, lg, idx, 1, 3, m[7], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 16, 0, 1, lg, m[7], 1, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 0, 1, lg, m[7], 1, 3, 1, qcf0);	//反向扫描
#endif
//9-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[4] * h[5];
	max = quanEVEN_DC(a + 16, 1, 0, lg, idx, 1, 3, m[8], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 16, 1, 0, lg, m[8], 1, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 1, 0, lg, m[8], 1, 3, 1, qcf0);	//反向扫描
#endif
//10-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = w[5] * h[5];
	max = quanEVEN_DC(a + 16, 1, 1, lg, idx, 1, 3, m[9], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 16, 1, 1, lg, m[9], 1, 3, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+16, 1, 1, lg, m[9], 1, 3, 1, qcf0);	//反向扫描
#endif

	/*第二次harr变换后的高频部分*/
//1-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[0] * h[0]) << 1;//第三次PTV
	max = quanEVEN_DC(a + 8, 0, 0, lg, idx, 3, 2, m[10], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 0, 0, lg, m[10], 3, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 0, 0, lg, m[10], 3, 2, 1, qcf0);	//反向扫描
#endif
//2-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[1] * h[0]) << 1;
	max = quanEVEN_DC(a + 8, 0, 4, lg, idx, 3, 2, m[11], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 0, 4, lg, m[11], 3, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 0, 4, lg, m[11], 3, 2, 1, qcf0);	//反向扫描
#endif
//3-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[0] * h[1]) << 1;
	max = quanEVEN_DC(a + 8, 4, 0, lg, idx, 3, 2, m[12], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 4, 0, lg, m[12], 3, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 4, 0, lg, m[12], 3, 2, 1, qcf0);	//反向扫描
#endif
//4-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[1] * h[1]) << 1;
	max = quanEVEN_DC(a + 8, 4, 4, lg, idx, 3, 2, m[13], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 4, 4, lg, m[13], 3, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 4, 4, lg, m[13], 3, 2, 1, qcf0);	//反向扫描
#endif
//5-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[3] * h[2]) << 1;//第二次PTV
	max = quanEVEN_DC(a + 8, 0, 2, lg, idx, 2, 2, m[14], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 0, 2, lg, m[14], 2, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 0, 2, lg, m[14], 2, 2, 1, qcf0);	//反向扫描
#endif
//6-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[2] * h[3]) << 1;
	max = quanEVEN_DC(a + 8, 2, 0, lg, idx, 2, 2, m[15], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 2, 0, lg, m[15], 2, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 2, 0, lg, m[15], 2, 2, 1, qcf0);	//反向扫描
#endif
//7-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[3] * h[3]) << 1;
	max = quanEVEN_DC(a + 8, 2, 2, lg, idx, 2, 2, m[16], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 2, 2, lg, m[16], 2, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 2, 2, lg, m[16], 2, 2, 1, qcf0);	//反向扫描
#endif
//8-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[5] * h[4]) << 1;//第一次PTV
	max = quanEVEN_DC(a + 8, 0, 1, lg, idx, 1, 2, m[17], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 0, 1, lg, m[17], 1, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 0, 1, lg, m[17], 1, 2, 1, qcf0);	//反向扫描
#endif
//9-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[4] * h[5]) << 1;
	max = quanEVEN_DC(a + 8, 1, 0, lg, idx, 1, 2, m[18], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 1, 0, lg, m[18], 1, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 1, 0, lg, m[18], 1, 2, 1, qcf0);	//反向扫描
#endif
//10-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[5] * h[5]) << 1;
	max = quanEVEN_DC(a + 8, 1, 1, lg, idx, 1, 2, m[19], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 8, 1, 1, lg, m[19], 1, 2, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+8, 1, 1, lg, m[19], 1, 2, 1, qcf0);	//反向扫描
#endif

	/*第一次harr变换后的高频部分*/
//1-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[2] * h[2]) << 2;	//第二次PTV
	max = quanEVEN_DC(a + 4, 0, 0, lg, idx, 2, 1, m[20], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 4, 0, 0, lg, m[20], 2, 1, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+4, 0, 0, lg, m[20], 2, 1, 1, qcf0);	//反向扫描
#endif
//2-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[3] * h[2]) << 2;
	max = quanEVEN_DC(a + 4, 0, 2, lg, idx, 2, 1, m[21], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 4, 0, 2, lg, m[21], 2, 1, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+4, 0, 2, lg, m[21], 2, 1, 1, qcf0);	//反向扫描
#endif
//3-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[2] * h[3]) << 2;
	max = quanEVEN_DC(a + 4, 2, 0, lg, idx, 2, 1, m[22], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 4, 2, 0, lg, m[22], 2, 1, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+4, 2, 0, lg, m[22], 2, 1, 1, qcf0);	//反向扫描
#endif

//4-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[3] * h[3]) << 2;
	max = quanEVEN_DC(a + 4, 2, 2, lg, idx, 2, 1, m[23], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 4, 2, 2, lg, m[23], 2, 1, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+4, 2, 2, lg, m[23], 2, 1, 1, qcf0);	//反向扫描
#endif

//5-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[5] * h[4]) << 2;	//第一次PTV
	max = quanEVEN_DC(a + 4, 0, 1, lg, idx, 1, 1, m[24], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 4, 0, 1, lg, m[24], 1, 1, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+4, 0, 1, lg, m[24], 1, 1, 1, qcf0);	//反向扫描
#endif

//6-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[4] * h[5]) << 2;
	max = quanEVEN_DC(a + 4, 1, 0, lg, idx, 1, 1, m[25], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 4, 1, 0, lg, m[25], 1, 1, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+4, 1, 0, lg, m[25], 1, 1, 1, qcf0);	//反向扫描
#endif
//7-------------------------------------------------------------------
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[5] * h[5]) << 2;
	max = quanEVEN_DC(a + 4, 1, 1, lg, idx, 1, 1, m[26], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);	//反量化
	scan_inv(a + 4, 1, 1, lg, m[26], 1, 1, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);	//反量化
	scan_inv(a+4, 1, 1, lg, m[26], 1, 1, 1, qcf0);	//反向扫描
#endif

	return;
}

void quant_subDC_noharr(int a, int *w, int *h, int **m, float delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int max, lg, idx = 0, qctr, maxcf0 = 0;
	extern union Fabs f1;
	float delta0;
	Uint8_Dat sign;
	delta0 = 1.0 / delta;

	DE_S_SUB desub;
	float *qcf0;
	int lenbits;
	int qnt = 0;
	/*PTV第二层*/
	/*1------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[2] * h[2]) << 3;
	max = quanEVEN_DC(a, 0, 0, lg, idx, 2, 0, m[27], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, 0, 0, lg, m[27], 2, 0, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 0, lg, m[27], 2, 0, 1, qcf0);	//反向扫描
#endif
	/*2------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[3] * h[2]) << 3;
	max = quanEVEN_DC(a, 0, 2, lg, idx, 2, 0, m[28], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, 0, 2, lg, m[28], 2, 0, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 2, lg, m[28], 2, 0, 1, qcf0);	//反向扫描
#endif
	/*3----------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[2] * h[3]) << 3;
	max = quanEVEN_DC(a, 2, 0, lg, idx, 2, 0, m[29], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, 2, 0, lg, m[29], 2, 0, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 2, 0, lg, m[29], 2, 0, 1, qcf0);	//反向扫描
#endif
	/*4--------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[3] * h[3]) << 3;
	max = quanEVEN_DC(a, 2, 2, lg, idx, 2, 0, m[30], delta*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, 2, 2, lg, m[30], 2, 0, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 2, 2, lg, m[30], 2, 0, 1, qcf0);	//反向扫描
#endif

	/*PTV第一层*/
	/*1----------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[5] * h[4]) << 3;
	max = quanEVEN_DC(a, 0, 1, lg, idx, 1, 0, m[31], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, 0, 1, lg, m[31], 1, 0, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 1, lg, m[31], 1, 0, 1, qcf0);	//反向扫描
#endif
	/*2-------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[4] * h[5]) << 3;
	max = quanEVEN_DC(a, 1, 0, lg, idx, 1, 0, m[32], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, 1, 0, lg, m[32], 1, 0, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 1, 0, lg, m[32], 1, 0, 1, qcf0);	//反向扫描
#endif
	/*3-----------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	lg = (w[5] * h[5]) << 3;
	max = quanEVEN_DC(a, 1, 1, lg, idx, 1, 0, m[33], delta0*10.0);
	qctr = quant_sub0(sn[0], max, nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, 1, 1, lg, m[33], 1, 0, 1, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 1, 1, lg, m[33], 1, 0, 1, qcf0);	//反向扫描
#endif

	return;
}

void quant_sub7_harr(int a, int *w, int *h, int **m, float delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	extern union Fabs f1;
	int *max;
	int idx = 0, lg, qctr, maxcf0 = 0;;
	Uint8_Dat sign;
	
	float delta0;
	max = (int*)calloc(7, sizeof(int));
	if (!max){
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = (float)1.0 / delta;
	DE_S_SUB desub;
	int lenbits;
	int qnt = 0;
	float *qcf0;
	/*第三次harr变换后的低频部分*/
	for (int i = 1; i < 8; i++)
		memset(sn[i], 0, len*sizeof(uchar));
	lg = w[6] * h[6];
	quanEVEN7(max, a, lg, idx, 3, m[34], delta0 * 10);

	/*1------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[0], nc[0], lg, idx, delta0, f1.qcf[1]); max[0] = 0;
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a, 0, 1, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 1, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*2---------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[1], nc[1], lg, idx, delta0, f1.qcf[2]); max[1] = 0;
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a, 0, 2, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 2, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*3-------------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[2], nc[2], lg, idx, delta0, f1.qcf[3]); max[2] = 0;
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a, 0, 3, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 3, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*4--------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[3], nc[3], lg, idx, delta0, f1.qcf[4]); max[3] = 0;
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a, 0, 4, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 4, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*5-----------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[4], nc[4], lg, idx, delta0, f1.qcf[5]); max[4] = 0;
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a, 0, 5, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 5, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*6-----------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[5], nc[5], lg, idx, delta0, f1.qcf[6]); max[5] = 0;
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a, 0, 6, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 6, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif

	/*7------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[6], nc[6], lg, idx, delta0, f1.qcf[7]); max[6] = 0;
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0); 
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a, 0, 7, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 7, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif

	/*第三次harr变换后的高频部分*/
	lg = w[6] * h[6];
	quanEVEN7(max, a + 16, lg, idx, 3, m[34], delta0 * 10);
	/*1--------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[0], nc[0], lg, idx, delta0, f1.qcf[1]); max[0] = 0;
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a + 16, 0, 1, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a+16, 0, 1, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*2---------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[1], nc[1], lg, idx, delta0, f1.qcf[2]); max[1] = 0;
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a + 16, 0, 2, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, 0, 2, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*3---------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[2], nc[2], lg, idx, delta0, f1.qcf[3]); max[2] = 0;
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a + 16, 0, 3, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, 0, 3, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*4----------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[3], nc[3], lg, idx, delta0, f1.qcf[4]); max[3] = 0;
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a + 16, 0, 4, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, 0, 4, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*5------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[4], nc[4], lg, idx, delta0, f1.qcf[5]); max[4] = 0;
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a + 16, 0, 5, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, 0, 5, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*6-------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[5], nc[5], lg, idx, delta0, f1.qcf[6]); max[5] = 0;
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a + 16, 0, 6, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, 0, 6, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*7---------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[6], nc[6], lg, idx, delta0, f1.qcf[7]); max[6] = 0;
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a + 16, 0, 7, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, 0, 7, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif

	/*第二次harr变换后的高频部分*/
	lg = (w[6] * h[6]) << 1;
	quanEVEN7(max, a + 8, lg, idx, 2, m[35], delta0 * 10);
	/*1-------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[0], nc[0], lg, idx, delta0, f1.qcf[1]); max[0] = 0;
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a + 8, 0, 1, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, 0, 1, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*2-------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[1], nc[1], lg, idx, delta0, f1.qcf[2]); max[1] = 0;
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a + 8, 0, 2, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, 0, 2, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*3----------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[2], nc[2], lg, idx, delta0, f1.qcf[3]); max[2] = 0;
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a + 8, 0, 3, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, 0, 3, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*4-------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[3], nc[3], lg, idx, delta0, f1.qcf[4]); max[3] = 0;
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a + 8, 0, 4, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, 0, 4, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*5-------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[4], nc[4], lg, idx, delta0, f1.qcf[5]); max[4] = 0;
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a + 8, 0, 5, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, 0, 5, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*6---------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[5], nc[5], lg, idx, delta0, f1.qcf[6]); max[5] = 0;
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a + 8, 0, 6, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, 0, 6, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*7----------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[6], nc[6], lg, idx, delta0, f1.qcf[7]); max[6] = 0;
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a + 8, 0, 7, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, 0, 7, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif

	/*第一次harr变换后的高频部分*/
	lg = (w[6] * h[6]) << 2;
	quanEVEN7(max, a + 4, lg, idx, 1, m[36], delta0 * 10);
	/*1-------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[0], nc[0], lg, idx, delta0, f1.qcf[1]);
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a + 4, 0, 1, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, 0, 1, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*2---------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[1], nc[1], lg, idx, delta0, f1.qcf[2]);
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a + 4, 0, 2, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, 0, 2, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*3----------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[2], nc[2], lg, idx, delta0, f1.qcf[3]);
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a + 4, 0, 3, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, 0, 3, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*4--------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[3], nc[3], lg, idx, delta0, f1.qcf[4]);
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a + 4, 0, 4, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, 0, 4, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*5----------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[4], nc[4], lg, idx, delta0, f1.qcf[5]);
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a + 4, 0, 5, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, 0, 5, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*6-------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[5], nc[5], lg, idx, delta0, f1.qcf[6]);
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a + 4, 0, 6, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, 0, 6, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*7-------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[6], nc[6], lg, idx, delta0, f1.qcf[7]);
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a + 4, 0, 7, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, 0, 7, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif

	free(max);
	return;
}

void quant_sub7_noharr(int a, int *w, int *h, int **m, float delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int idx = 0, lg, qctr, maxcf0 = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	float delta0;
	max = (int*)calloc(7, sizeof(int));
	if (!max){
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = (float)1.0 / delta;
	DE_S_SUB desub;
	int lenbits;
	int qnt = 0;
	float *qcf0;
	/*DC后七个系数*/
	lg = (w[6] * h[6]) << 3;
	quanEVEN7(max, a, lg, idx, 0, m[37], delta0*10.0);
	/*1-------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[0], nc[0], lg, idx, delta0, f1.qcf[1]);
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a, 0, 1, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 1, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*2------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[1], nc[1], lg, idx, delta0, f1.qcf[2]);
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a, 0, 2, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 2, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*3-----------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[2], nc[2], lg, idx, delta0, f1.qcf[3]);
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a, 0, 3, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 3, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*4-------------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[3], nc[3], lg, idx, delta0, f1.qcf[4]);
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a, 0, 4, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 4, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*5--------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[4], nc[4], lg, idx, delta0, f1.qcf[5]);
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a, 0, 5, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 5, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*6------------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[5], nc[5], lg, idx, delta0, f1.qcf[6]);
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a, 0, 6, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 6, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*7--------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[6], nc[6], lg, idx, delta0, f1.qcf[7]);
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a, 0, 7, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, 0, 7, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif

	free(max);
	return;
}


 void quant_sub8_harr(int a, int b, int *w, int *h, int **m, float delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	float delta0;
	max = (int*)calloc(8, sizeof(int));
	if (!max){
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = (float)1.0 / delta;
	DE_S_SUB desub;
	int lenbits;
	int qnt = 0;
	float *qcf0;
	/*第三次harr变换后的低频部分*/
	lg = w[6] * h[6];
	quanEVEN8(max, a, b, lg, idx, 3, m[34], delta0*10.0);

	/*1------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[0], max[0], nc[0], lg, idx, delta0, f1.qcf[0]); max[0] = 0;
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, b, 0, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 0, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*2--------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[1], nc[1], lg, idx, delta0, f1.qcf[1]); max[1] = 0;
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a, b, 1, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 1, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*3---------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[2], nc[2], lg, idx, delta0, f1.qcf[2]); max[2] = 0;
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a, b, 2, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 2, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*4---------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[3], nc[3], lg, idx, delta0, f1.qcf[3]); max[3] = 0;
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a, b, 3, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 3, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*5-----------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[4], nc[4], lg, idx, delta0, f1.qcf[4]); max[4] = 0;
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a, b, 4, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 4, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*6-----------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[5], nc[5], lg, idx, delta0, f1.qcf[5]); max[5] = 0;
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a, b, 5, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 5, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*7----------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[6], nc[6], lg, idx, delta0, f1.qcf[6]); max[6] = 0;
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a, b, 6, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 6, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*8------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[7], nc[7], lg, idx, delta0, f1.qcf[7]); max[7] = 0;
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a, b, 7, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 7, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif

	/*第三次harr变换后的高频部分*/
	lg = w[6] * h[6];
	quanEVEN8(max, a + 16, b, lg, idx, 3, m[34], delta0*10.0);
	/*1--------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[0], max[0], nc[0], lg, idx, delta0, f1.qcf[0]); max[0] = 0;
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a + 16, b, 0, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, b, 0, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*2---------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[1], nc[1], lg, idx, delta0, f1.qcf[1]); max[1] = 0;
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a + 16, b, 1, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, b, 1, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*3--------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[2], nc[2], lg, idx, delta0, f1.qcf[2]); max[2] = 0;
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a + 16, b, 2, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, b, 2, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*4--------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[3], nc[3], lg, idx, delta0, f1.qcf[3]); max[3] = 0;
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a + 16, b, 3, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, b, 3, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*5--------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[4], nc[4], lg, idx, delta0, f1.qcf[4]); max[4] = 0;
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a + 16, b, 4, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, b, 4, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*6-------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[5], nc[5], lg, idx, delta0, f1.qcf[5]); max[5] = 0;
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a + 16, b, 5, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, b, 5, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*7--------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[6], nc[6], lg, idx, delta0, f1.qcf[6]); max[6] = 0;
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a + 16, b, 6, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, b, 6, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif
	/*8------------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[7], nc[7], lg, idx, delta0, f1.qcf[7]); max[7] = 0;
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a + 16, b, 7, lg, m[34], 3, 3, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 16, b, 7, lg, m[34], 3, 3, 0, qcf0);	//反向扫描
#endif

	/*第二次harr变换后的高频部分*/
	lg = (w[6] * h[6]) << 1;
	quanEVEN8(max, a + 8, b, lg, idx, 2, m[35], delta0*10.0);
	/*1-----------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[0], max[0], nc[0], lg, idx, delta0, f1.qcf[0]); max[0] = 0;
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a + 8, b, 0, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, b, 0, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*2-----------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[1], nc[1], lg, idx, delta0, f1.qcf[1]); max[1] = 0;
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a + 8, b, 1, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, b, 1, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*3-------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[2], nc[2], lg, idx, delta0, f1.qcf[2]); max[2] = 0;
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a + 8, b, 2, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, b, 2, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*4--------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[3], nc[3], lg, idx, delta0, f1.qcf[3]); max[3] = 0;
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a + 8, b, 3, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, b, 3, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*5----------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[4], nc[4], lg, idx, delta0, f1.qcf[4]); max[4] = 0;
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a + 8, b, 4, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, b, 4, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*6---------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[5], nc[5], lg, idx, delta0, f1.qcf[5]); max[5] = 0;
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a + 8, b, 5, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, b, 5, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*7------------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[6], nc[6], lg, idx, delta0, f1.qcf[6]); max[6] = 0;
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a + 8, b, 6, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, b, 6, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif
	/*8----------------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[7], nc[7], lg, idx, delta0, f1.qcf[7]); max[7] = 0;
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a + 8, b, 7, lg, m[35], 3, 2, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 8, b, 7, lg, m[35], 3, 2, 0, qcf0);	//反向扫描
#endif

	/*第一次harr变换后的高频部分*/
	lg = (w[6] * h[6]) << 2;
	quanEVEN8(max, a + 4, b, lg, idx, 1, m[36], delta0*10.0);
	/*1----------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[0], max[0], nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a + 4, b, 0, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, b, 0, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*2--------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[1], nc[1], lg, idx, delta0, f1.qcf[1]);
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a + 4, b, 1, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, b, 1, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*3--------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[2], nc[2], lg, idx, delta0, f1.qcf[2]);
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a + 4, b, 2, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, b, 2, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*4------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[3], nc[3], lg, idx, delta0, f1.qcf[3]);
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a + 4, b, 3, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, b, 3, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*5-------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[4], nc[4], lg, idx, delta0, f1.qcf[4]);
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a + 4, b, 4, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, b, 4, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*6-----------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[5], nc[5], lg, idx, delta0, f1.qcf[5]);
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a + 4, b, 5, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, b, 5, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*7-----------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[6], nc[6], lg, idx, delta0, f1.qcf[6]);
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a + 4, b, 6, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, b, 6, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif
	/*8------------------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[7], nc[7], lg, idx, delta0, f1.qcf[7]);
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a + 4, b, 7, lg, m[36], 3, 1, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a + 4, b, 7, lg, m[36], 3, 1, 0, qcf0);	//反向扫描
#endif

	free(max);

	return;
}

void quant_sub8_noharr(int a, int b, int *w, int *h, int **m, float delta)
{
	extern uchar **sn;
	extern uchar *snbin;
	extern uchar *bin;
	extern int ptr;
	extern int **nc;
	extern int len;
	int *max;
	int lg, idx = 0, qctr, maxcf0 = 0;
	extern union Fabs f1;
	Uint8_Dat sign;
	float delta0;
	float *qcf0;
	max = (int*)calloc(8, sizeof(int));
	if (!max){
		printf("创建max数组失败！\n");
		exit(1);
	}
	delta0 = (float)1.0 / delta;
	DE_S_SUB desub;
	int lenbits;
	int qnt = 0;

	lg = (w[6] * h[6]) << 3;
	quanEVEN8(max, a, b, lg, idx, 0, m[37], delta0*10.0);
	/*1-------------------------------------------------------*/
	memset(bin, 0, len*sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[0], max[0], nc[0], lg, idx, delta0, f1.qcf[0]);
#if ENCODE
	sign.len = handleSn(sn[0], idx, lg, qctr, &maxcf0, f1.qcf[0]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[0], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[0], sn[0], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[0], delta, lg);
	scan_inv(a, b, 0, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 0, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*2-----------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[1], max[1], nc[1], lg, idx, delta0, f1.qcf[1]);
#if ENCODE
	sign.len = handleSn(sn[1], idx, lg, qctr, &maxcf0, f1.qcf[1]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[1], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[1], sn[1], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[1], delta, lg);
	scan_inv(a, b, 1, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 1, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*3------------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[2], max[2], nc[2], lg, idx, delta0, f1.qcf[2]);
#if ENCODE
	sign.len = handleSn(sn[2], idx, lg, qctr, &maxcf0, f1.qcf[2]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[2], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[2], sn[2], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[2], delta, lg);
	scan_inv(a, b, 2, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 2, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*4---------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[3], max[3], nc[3], lg, idx, delta0, f1.qcf[3]);
#if ENCODE
	sign.len = handleSn(sn[3], idx, lg, qctr, &maxcf0, f1.qcf[3]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[3], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[3], sn[3], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[3], delta, lg);
	scan_inv(a, b, 3, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 3, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*5--------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[4], max[4], nc[4], lg, idx, delta0, f1.qcf[4]);
#if ENCODE
	sign.len = handleSn(sn[4], idx, lg, qctr, &maxcf0, f1.qcf[4]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[4], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[4], sn[4], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[4], delta, lg);
	scan_inv(a, b, 4, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 4, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*6-------------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[5], max[5], nc[5], lg, idx, delta0, f1.qcf[5]);
#if ENCODE
	sign.len = handleSn(sn[5], idx, lg, qctr, &maxcf0, f1.qcf[5]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[5], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[5], sn[5], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[5], delta, lg);
	scan_inv(a, b, 5, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 5, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*7----------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[6], max[6], nc[6], lg, idx, delta0, f1.qcf[6]);
#if ENCODE
	sign.len = handleSn(sn[6], idx, lg, qctr, &maxcf0, f1.qcf[6]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[6], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[6], sn[6], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[6], delta, lg);
	scan_inv(a, b, 6, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 6, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif
	/*8----------------------------------------------------------------------*/
	memset(bin, 0, len * sizeof(uchar));
	idx = 0; ptr = 0;
	qctr = quant_sub0(sn[7], max[7], nc[7], lg, idx, delta0, f1.qcf[7]);
#if ENCODE
	sign.len = handleSn(sn[7], idx, lg, qctr, &maxcf0, f1.qcf[7]);
	sign.dat = snbin;
	en_sub3d_sub2(f1.qcf[7], &sign, lg, maxcf0);
#endif
#if DECODE
	lenbits = ptr;
	ptr = 0;
	qctr = find_qctr(bin, &qnt);
	if (qnt){
		desub = de_sub3d_sub2(bin, lg, qctr, lenbits);	//解码
		free(desub.sn.dat); free(desub.cfw.dat);
	}
	else{
		printf("子带全为零\n");
	}
#endif
#if NOENCODE
	addSn(f1.qcf[7], sn[7], lg);
	ptr = 0;
	qcf0 = rst_sub(bin, f1.qcf[7], delta, lg);
	scan_inv(a, b, 7, lg, m[37], 3, 0, 0, qcf0);
#endif
#if INVQUANT
	ptr = 0;
	qcf0 = rst_sub(bin, desub.cf.dat, delta, lg);
	scan_inv(a, b, 7, lg, m[37], 3, 0, 0, qcf0);	//反向扫描
#endif

	free(max);

	return;
}



/*输入为W,H; W为宽，H为高，例如1024*1024*96的数据，W=1024，H=1024*
*输出为w,h,m; w和h分别为PTV每个小块的宽和高，m为希尔伯特扫描路径*
*w,h为大小为7的一维数组，m是38行的二维数组，对应38个希尔伯特路径*/
int GetM(int W, int H, int *w, int *h, int **m)
{
	int h1h, w1h, h2h, w2h, h3h, w3h;
	int h3l, w3l, h2l, w2l, h1l, w1l;
	W >>= 3; H >>= 3;
	/*求PTV后每个小块的大小*/
	if (!(W % 2)){
		w3h = W >> 1; w3l = W >> 1;
	}
	else{
		w3h = (W + 1) >> 1; w3l = w3h - 1;
	}
	if (!(w3h % 2)){
		w2h = w3h >> 1; w2l = w3h >> 1;
	}
	else{
		w2h = (w3h + 1) >> 1; w2l = w2h - 1;
	}
	if (!(w2h % 2)){
		w1h = w2h >> 1; w1l = w1h;
	}
	else{
		w1h = (w2h + 1) >> 1; w1l = w1h - 1;
	}
	if (!(H % 2)){
		h3h = H >> 1; h3l = H >> 1;
	}
	else{
		h3h = (H + 1) >> 1; h3l = h3h - 1;
	}
	if (!(h3h % 2)){
		h2h = h3h >> 1; h2l = h2h;
	}
	else{
		h2h = (h3h + 1) >> 1; h2l = h2h - 1;
	}
	if (!(h2h % 2)){
		h1h = h2h >> 1; h1l = h1h;
	}
	else{
		h1h = (h2h + 1) >> 1; h1l = h1h - 1;
	}
	w[0] = w1h; w[1] = w1l;
	w[2] = w2h; w[3] = w2l;
	w[4] = w3h; w[5] = w3l;
	w[6] = W;
	h[0] = h1h; h[1] = h1l;
	h[2] = h2h; h[3] = h2l;
	h[4] = h3h; h[5] = h3l;
	h[6] = H;
	
	/*求路径M*/
	m[0] = modify_Hilbert2D(w1h, h1h);
	m[1] = modify_Hilbert2D(w1l, h1h);
	m[2] = modify_Hilbert2D(w1h, h1l);
	m[3] = modify_Hilbert2D(w1l, h1l);

	m[4] = modify_Hilbert2D(w2l, h2h);
	m[5] = modify_Hilbert2D(w2h, h2l);
	m[6] = modify_Hilbert2D(w2l, h2l);

	m[7] = modify_Hilbert2D(w3l, h3h);
	m[8] = modify_Hilbert2D(w3h, h3l);
	m[9] = modify_Hilbert2D(w3l, h3l);
	/**/
	m[10] = modify_Hilbert3D(w1h, h1h, 1);
	m[11] = modify_Hilbert3D(w1l, h1h, 1);
	m[12] = modify_Hilbert3D(w1h, h1l, 1);
	m[13] = modify_Hilbert3D(w1l, h1l, 1);

	m[14] = modify_Hilbert3D(w2l, h2h, 1);
	m[15] = modify_Hilbert3D(w2h, h2l, 1);
	m[16] = modify_Hilbert3D(w2l, h2l, 1);

	m[17] = modify_Hilbert3D(w3l, h3h, 1);
	m[18] = modify_Hilbert3D(w3h, h3l, 1);
	m[19] = modify_Hilbert3D(w3l, h3l, 1);
	/**/
	m[20] = modify_Hilbert3D(w2h, h2h, 2);//32*32*4
	m[21] = modify_Hilbert3D(w2l, h2h, 2);
	m[22] = modify_Hilbert3D(w2h, h2l, 2);
	m[23] = modify_Hilbert3D(w2l, h2l, 2);

	m[24] = modify_Hilbert3D(w3l, h3h, 2);//64*64*4
	m[25] = modify_Hilbert3D(w3h, h3l, 2);
	m[26] = modify_Hilbert3D(w3l, h3l, 2);
	/**/
	m[27] = modify_Hilbert3D(w2h, h2h, 3);//32*32*8
	m[28] = modify_Hilbert3D(w2l, h2h, 3);
	m[29] = modify_Hilbert3D(w2h, h2l, 3);
	m[30] = modify_Hilbert3D(w2l, h2l, 3);

	m[31] = modify_Hilbert3D(w3l, h3h, 3);//64*64*8
	m[32] = modify_Hilbert3D(w3h, h3l, 3);
	m[33] = modify_Hilbert3D(w3l, h3l, 3);

	m[34] = modify_Hilbert2D(W, H);//128*128
	m[35] = modify_Hilbert3D(W, H, 1);//128*128*2
	m[36] = modify_Hilbert3D(W, H, 2);//128*128*4
	m[37] = modify_Hilbert3D(W, H, 3);//128*128*8
	return 0;
}

int handleSn(uchar *sgn, int idx, int lg, int qctr, int *maxcf, float *fabsA)
{
	extern unsigned char *snbin;
	//unsigned char *snbin = (unsigned char*)malloc((lg >> 3)*sizeof(unsigned char));
	unsigned char temp = 0;
	int cnt = 7, num = 0, max = 0;
	int ia = 0;
	int len_snbin;
	if (qctr){ //quantTHD量化的
		for (int i = idx; i < idx + lg; i++){
			if (fabsA[i]>max)
				max = fabsA[i];
			if (sgn[i]){
				sgn[i]--;
				temp |= (sgn[i] << cnt);
				cnt--;
				ia++;
			}
			if (cnt == -1){
				cnt = 7;
				snbin[num] = temp;
				num++;
				temp = 0;
			}
		}
	}
	else{ //quantEVEN量化，这里量化后是没有0值的。因为是向上取整
		for (int i = idx; i < idx + lg; i++){
			fabsA[i]--;
			if (fabsA[i]>max)
				max = fabsA[i];
			if (!sgn[i]){
				temp <<= 1;
				cnt--;
			}
			else{
				sgn[i]--;
				temp |= (sgn[i] << cnt);
				cnt--;
			}
			if (cnt == -1){
				cnt = 7;
				snbin[num] = temp;
				num++;
				temp = 0;
			}
		}
	}
	if (cnt == 7){
		len_snbin = num * 8;
	}
	else{
		snbin[num] = temp; num++;
		len_snbin = num * 8 - cnt - 1;
	}
	*maxcf = max;
	return len_snbin; //此处的长度是snbin的总位数，不是数组的长度。
}