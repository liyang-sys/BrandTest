#include "quant_inv.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "all.h"

float* rst_sub(uchar* bin, int* qcf, float delta, int lg)
//float* rst_sub(uchar* bin, float* cf, float delta, int lg)
{
	extern int len;
	extern int ptr;
	extern float ***rePTVData;
	int max=0, temp;
	int *nc = (int*)calloc(len, sizeof(int));
	char *sgn = (char*)calloc(lg, sizeof(char));
	float *qf = NULL, T;
	uchar temp1 = bin[0] & 0x80; temp1 >>= 7;
	uchar temp2 = bin[0] & 0x40; temp2 >>= 6;
	uchar temp3 = bin[0] & 0x20; temp3 >>= 5;

	//int *qcf = (int*)malloc(lg*sizeof(int));	//临时添加
	//for (int i = 0; i < lg; i++)
	//	qcf[i] = (int)cf[i];

	if (temp1){
		if (temp2){	//trim
			T = 0.55 * delta;
			for (int i = 0; i < lg; i++){
				if (qcf[i] & 0x80000000){
					sgn[i]--;
					qcf[i] *= -1;
				}
				else
					sgn[i]++;
				temp = qcf[i] + 1;
				if (temp>max)
					max = temp;
				nc[temp]++;
			}
			qf = rstTHDctr(qcf, sgn, T, delta, nc, max, lg);
		}
		else{
			if (temp3){	//quant center
				for (int i = 0; i < lg; i++){
					if (qcf[i] & 0x80000000){
						sgn[i]--;
						qcf[i] *= -1;
					}
					else
						sgn[i]++;
					temp = qcf[i] + 1;
					if (temp>max)
						max = temp;
					nc[temp]++;
				}
				qf = rstTHDctr(qcf, sgn, 0.52*delta, delta, nc, max, lg);
			}
			else{	//quant even
				ptr += 3;
				DES de = deSFcode(bin, 64);
				for (int i = 0; i < lg; i++){
					if (qcf[i] & 0x80000000){
						sgn[i]--;
						qcf[i] *= -1;
					}
					else
						sgn[i]++;
					temp = qcf[i];
					if (temp>max)
						max = temp;
					nc[temp]++;
				}
				qf = rstEVENctr(qcf, sgn, delta, nc, max, (float)de.sym-1, lg);
			}
		}
	}
	free(nc);
	free(qcf); free(sgn);
	return qf;
}

float* rstTHDctr(int* qcf, char *sgn, float T, float delta, int* nc, int maxcf, int lg)
{
	float a3 = 1.0 / 3.0, b3 = 2.0 / 3.0, na, nb;
	int maxval = maxcf - 1;
	int valT = 1, val, idx=1;
	int minnc = nc[1];
	float* qf = (float*)malloc(lg*sizeof(float));
	for (int i = 2; i < maxcf + 1; i++){
		if (nc[i] < minnc){
			minnc = nc[i];
			idx = i;
		}
	}
	
	for (int i = 1; i < idx - 1; i++){
		na = nc[i + 1]; nb = nc[i + 2];
		if (nb>na){
			valT = i;
			break;
		}
		else
			valT = i + 1;
	}
	
	for (int i = 0; i < lg; i++){
		val = qcf[i];
		if (!val){
			qf[i] = 0;
		}
		else if (val < valT){
			na = nc[val + 1]; nb = nc[val + 2];
			float ctr_geo = (a3*na + b3*nb) / (na + nb)*delta;
			qf[i] = (T + ((float)val - 1)*delta + ctr_geo);
		}
		else{
			if (maxval>1){
				qf[i] = (T + (val - 0.5)*delta);
			}
			else{
				qf[i] = (T + 0.37*delta);
			}
		}
		qf[i] *= sgn[i];
	}

	return qf;
}

float* rstEVENctr(int* qcf, char *sgn, float delta, int* nc, int maxcf, float ctr1, int lg)
{
	float na, nb, a3 = 1.0 / 3.0, b3 = 2.0 / 3.0;
	int maxval = maxcf, minnc = nc[1], val, idx = 1;
	float* qf = (float*)malloc(lg*sizeof(float));

	for (int i = 2; i < maxcf + 1; i++){
		if (nc[i] < minnc){
			minnc = nc[i];
			idx = i;
		}
	}
	ctr1 = (float)((float)ctr1*0.125 / 64.0 + 0.375)*delta;
	for (int i = 0; i < lg; i++){
		val = qcf[i];
		qf[i] = ((float)val - 1.0)*delta;
		if (val == 1){
			qf[i] += ctr1;
		}
		else if (val < idx - 1){
			na = nc[val]; nb = nc[val + 1];
			float ctr_geo = (a3*na + b3*nb) / (na + nb)*delta;
			qf[i] += ctr_geo;
		}
		else{
			qf[i] += (0.5*delta);
		}
		qf[i] *= sgn[i];
	}


	return qf;
}

void quant_DCinv(char* bin, int a, int lg, int* m, float delta)
{
	union data {
		unsigned long long a;
		uchar b[8];
	} rem1, rem2;
	rem1.a = 0; rem2.a = 0;
	extern int ptr;
	extern float ***rePTVData;
	int nbits, maxqsb, minqsb, nbitsMax, Nsym;
	int b = 0, c = 0;
	int idx = ptr >> 3;
	rem1.b[5] = bin[idx];
	rem1.b[4] = bin[idx + 1];
	rem1.b[3] = bin[idx + 2];
	rem1.b[2] = bin[idx + 3];
	rem1.b[1] = bin[idx + 4];
	rem1.b[0] = bin[idx + 5];
	/*解码maxqsb*/
	nbits = (int)(17.5 - log(delta) / log(2.0));
	rem1.a <<= nbits;
	rem2.b[1] = rem1.b[7];
	rem2.b[0] = rem1.b[6];
	maxqsb = (int)rem2.a;
	rem1.b[7] = 0;
	rem1.b[6] = 0;
	ptr += nbits;
	/*解码minqsb*/
	nbitsMax = (int)(log(maxqsb) / log(2)) + 1;
	rem1.a <<= nbitsMax;
	rem2.b[1] = rem1.b[7];
	rem2.b[0] = rem1.b[6];
	minqsb = (int)rem2.a;
	rem1.b[7] = 0;
	rem1.b[6] = 0;
	ptr += nbitsMax;

	Nsym = maxqsb - minqsb + 1;

	/*解码bin*/
	int j = 0;
	int offset0 = 8;
	int offset1 = 32;
	for (int i = 0; i < lg; i++){
		DES de = deSFcode(bin, Nsym);
		rePTVData[a][b][c] = (float)((de.sym - 1 + minqsb) * delta);
		switch (m[j])
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
		}j++;
	}

	return;
}

void scan_inv(int a, int b, int c, int lg, int *m, int ptv, int harr, int dc, float *qcf)
{
	extern float ***rePTVData;
	extern float ***reVideoData;
	float ***data;
	int offset0, offset1, i = 0;
	if (!qcf){
		printf("qcf为空，或qcf全为零\n");
		return;
	}
	if (dc){
		data = rePTVData;
		offset0 = 1 << ptv;
	}
	else{
		data = reVideoData;
		offset0 = 8;
	}
	offset1 = 4 << harr;
	for (int j = 0; j < lg; j++){
		data[a][b][c] = qcf[j];
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
	free(qcf);
	
	return;
}

int find_qctr(uchar *bin, int *qnt)
{
	union data {
		unsigned long long a;
		uchar b[8];
	} rem1;
	extern int ptr;
	int qctr;
	rem1.a = 0;
	uchar qnt0;
	int idx = ptr >> 3;
	rem1.b[5] = bin[idx];
	rem1.b[4] = bin[idx + 1];
	rem1.b[3] = bin[idx + 2];
	rem1.b[2] = bin[idx + 3];
	rem1.b[1] = bin[idx + 4];
	rem1.b[0] = bin[idx + 5];

	rem1.a <<= 1;
	qnt0 = (uchar)rem1.b[6];
	rem1.b[6] = 0;

	*qnt = qnt0;
	if (qnt0){
		uchar trim;
		rem1.a <<= 1;
		trim = (uchar)rem1.b[6];
		rem1.b[6] = 0;
		if (trim){
			qctr = 1;
			ptr += 2;
		}
		else{
			rem1.a <<= 1;
			qctr = (int)rem1.b[6];
			rem1.b[6] = 0;
			ptr += 3;
			if (!qctr){
				deSFcode(bin, 64); //这里可以优化，这里解码的ctr1后面量化时需要，所以后面量化时又重复了这一步
			}
		}
	}
	else{
		qctr = 1;
	}

	return qctr;
}