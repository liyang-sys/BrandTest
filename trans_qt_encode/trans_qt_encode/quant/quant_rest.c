#include "quant_rest.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

union Fabs{
	float **qcf;
	unsigned int **temp;
};

int quanEVEN_DC(int a, int b, int c, int lg, int pos, int ptv, int harr, int *m, float delta)
{
	extern unsigned char **sn;
	extern int **nc;
	int i = 0, j, lth, qf, max=0, offset1, offset0;
	extern float ***PTVData;
	extern union Fabs f1;
	lth = pos + lg;
	offset0 = 1 << ptv;
	offset1 = 4 << harr;
	for (j = pos; j < lth; j++){
		f1.qcf[0][j] = PTVData[a][b][c];
		if (f1.qcf[0][j]){
			sn[0][j] = (((f1.temp[0][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[0][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[0][j] &= 0x7fffffff; //������λ���������ֵ
		qf = (int)(f1.qcf[0][j] * delta);
		if (qf > max){
			max = qf;
		}
		nc[0][qf]++;//ͳ��nc
		switch (m[i])
		{
		case 3:
			c+=offset0;
			break;
		case -3:
			c-=offset0;
			break;
		case 2:
			b+=offset0;
			break;
		case -2:
			b-=offset0;
			break;
		case 1:
			a+=offset1;
			break;
		case -1:
			a-=offset1;
			break;
		default:
			break;
		}i++;
	}
	return max;
}

int quanEVEN8(int *max, int a, int b, int lg, int pos, int harr, int *m, float delta0) //b,a�ֱ�Ϊ�����к���,����8���Ӵ�
{
	extern unsigned char **sn;
	int i=0, j, qf, c = 0, offset;
	extern int **nc;
	extern float ***VideoData;
	extern union Fabs f1;
	lg += pos;
	offset = 4 << harr;
	for (j = pos; j < lg; j++){
		/*��1���Ӵ�*/
		f1.qcf[0][j] = VideoData[a][b][c];
		if (f1.qcf[0][j]){
			sn[0][j] = (((f1.temp[0][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[0][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[0][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[0][j] * delta0);	//����
		if (qf > max[0]){//���������ϵ�������ֵ
			max[0] = qf;
		}
		nc[0][qf]++; //��nc

		/*��2���Ӵ�*/
		f1.qcf[1][j] = VideoData[a][b][c + 1];
		if (f1.qcf[1][j]){
			sn[1][j] = (((f1.temp[1][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[1][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[1][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[1][j] * delta0);	//����
		if (qf > max[1]){//���������ϵ�������ֵ
			max[1] = qf;
		}
		nc[1][qf]++; //��nc

		/*3*/
		f1.qcf[2][j] = VideoData[a][b][c + 2];
		if (f1.qcf[2][j]){
			sn[2][j] = (((f1.temp[2][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[2][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[2][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[2][j] * delta0);	//����
		if (qf > max[2]){//���������ϵ�������ֵ
			max[2] = qf;
		}
		nc[2][qf]++; //��nc

		/*4*/
		f1.qcf[3][j] = VideoData[a][b][c + 3];
		if (f1.qcf[3][j]){
			sn[3][j] = (((f1.temp[3][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[3][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[3][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[3][j] * delta0);	//����
		if (qf > max[3]){//���������ϵ�������ֵ
			max[3] = qf;
		}
		nc[3][qf]++; //��nc

		/*5*/
		f1.qcf[4][j] = VideoData[a][b][c + 4];
		if (f1.qcf[4][j]){
			sn[4][j] = (((f1.temp[4][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[4][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[4][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[4][j] * delta0);	//����
		if (qf > max[4]){//���������ϵ�������ֵ
			max[4] = qf;
		}
		nc[4][qf]++; //��nc

		/*6*/
		f1.qcf[5][j] = VideoData[a][b][c + 5];
		if (f1.qcf[5][j]){
			sn[5][j] = (((f1.temp[5][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[5][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[5][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[5][j] * delta0);	//����
		if (qf > max[5]){//���������ϵ�������ֵ
			max[5] = qf;
		}
		nc[5][qf]++; //��nc

		/*7*/
		f1.qcf[6][j] = VideoData[a][b][c + 6];
		if (f1.qcf[6][j]){
			sn[6][j] = (((f1.temp[6][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[6][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[6][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[6][j] * delta0);	//����
		if (qf > max[6]){//���������ϵ�������ֵ
			max[6] = qf;
		}
		nc[6][qf]++; //��nc

		/*8*/
		f1.qcf[7][j] = VideoData[a][b][c + 7];
		if (f1.qcf[7][j]){
			sn[7][j] = (((f1.temp[7][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[7][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[7][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[7][j] * delta0);	//����
		if (qf > max[7]){//���������ϵ�������ֵ
			max[7] = qf;
		}
		nc[7][qf]++; //��nc

		/*����sn*/
		switch (m[i])
		{
		case 1:
			a += offset;
			break;
		case -1:
			a -= offset;
			break;
		case 2:
			b += 8;
			break;
		case -2:
			b -= 8;
			break;
		case 3:
			c += 8;
			break;
		case -3:
			c -= 8;
			break;
		default:
			break;
		}i++;
	}
	return lg;
}

int quanEVEN7(int *max, int a, int lg, int pos, int harr, int *m, float delta0) //aΪ������
{
	extern unsigned char **sn;
	extern int **nc;
	int i = 0, j, qf, c = 1, b = 0, offset;
	extern float ***VideoData;
	extern union Fabs f1;
	lg += pos;
	offset = 4 << harr;

	for (j = pos; j < lg; j++){
		/*1*/
		f1.qcf[1][j] = VideoData[a][b][c];
		if (f1.qcf[1][j]){
			sn[1][j] = (((f1.temp[1][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[1][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[1][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[1][j] * delta0);	//����
		if (qf > max[0]){//���������ϵ�������ֵ
			max[0] = qf;
		}
		nc[0][qf]++; //��nc

		/*2*/
		f1.qcf[2][j] = VideoData[a][b][c + 1];
		if (f1.qcf[2][j]){
			sn[2][j] = (((f1.temp[2][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[2][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[2][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[2][j] * delta0);	//����
		if (qf > max[1]){//���������ϵ�������ֵ
			max[1] = qf;
		}
		nc[1][qf]++; //��nc

		/*3*/
		f1.qcf[3][j] = VideoData[a][b][c + 2];
		if (f1.qcf[3][j]){
			sn[3][j] = (((f1.temp[3][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[3][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[3][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[3][j] * delta0);	//����
		if (qf > max[2]){//���������ϵ�������ֵ
			max[2] = qf;
		}
		nc[2][qf]++; //��nc

		/*4*/
		f1.qcf[4][j] = VideoData[a][b][c + 3];
		if (f1.qcf[4][j]){
			sn[4][j] = (((f1.temp[4][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[4][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[4][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[4][j] * delta0);	//����
		if (qf > max[3]){//���������ϵ�������ֵ
			max[3] = qf;
		}
		nc[3][qf]++; //��nc

		/*5*/
		f1.qcf[5][j] = VideoData[a][b][c + 4];
		if (f1.qcf[5][j]){
			sn[5][j] = (((f1.temp[5][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[5][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[5][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[5][j] * delta0);	//����
		if (qf > max[4]){//���������ϵ�������ֵ
			max[4] = qf;
		}
		nc[4][qf]++; //��nc

		/*6*/
		f1.qcf[6][j] = VideoData[a][b][c + 5];
		if (f1.qcf[6][j]){
			sn[6][j] = (((f1.temp[6][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[6][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[6][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[6][j] * delta0);	//����
		if (qf > max[5]){//���������ϵ�������ֵ
			max[5] = qf;
		}
		nc[5][qf]++; //��nc

		/*7*/
		f1.qcf[7][j] = VideoData[a][b][c + 6];
		if (f1.qcf[7][j]){
			sn[7][j] = (((f1.temp[7][j] & 0x80000000) >> 31) & 0x01);//ȡ����
			sn[7][j]++; //1��ʾ����2��ʾ����0��ʾ��
		}
		f1.temp[7][j] &= 0x7fffffff;		//ͨ��������λ�����������ֵ
		qf = (int)(f1.qcf[7][j] * delta0);	//����
		if (qf > max[6]){//���������ϵ�������ֵ
			max[6] = qf;
		}
		nc[6][qf]++; //��nc

		/*����sn*/
		switch (m[i])
		{
		case 3:
			c += 8;
			break;
		case -3:
			c -= 8;
			break;
		case 2:
			b += 8;
			break;
		case -2:
			b -= 8;
			break;
		case 1:
			a += offset;
			break;
		case -1:
			a -= offset;
			break;
		default:
			break;
		}i++;
	}
	return lg;
}

int quanEVEN2(int len, int pos, float *absA, float delta)  //ֻ��������������nc��runs
{
	int i, num1 = 0;
	int ctr[2];
	float sum1 = 0.0;
	float temp;
	len += pos;
	for (i = pos; i < len; i++){
		temp = *(absA + i);
		*(absA + i) = (int)(temp * delta) + 1;
		if (*(absA + i) == 1){
			num1++;
			sum1 += temp;
		}
	}
	if (num1 > 0){
		ctr[0] = (int)(((sum1*delta) / num1) * 512 - 191.5);
	}
	else ctr[0] = 63;
	return ctr[0];
}

QT quanTHD(uchar *sgn, int len, int pos, float *fabsA, float delta)
{
	QT qtTHD;
	int i, j = 0, *temp;
	extern float *absA2;
	extern float *absA;
	qtTHD.runs = (int*)calloc(len, sizeof(int));
	if (!qtTHD.runs){
		printf("������̬����runsʧ�ܣ�\n");
		exit(1);
	}
	len += pos;
	for (i = pos; i < len; i++)
	{
		qtTHD.runs[j]++;		//��runs
		absA2[i] = fabsA[i] * fabsA[i];
		absA[i] = fabsA[i];
		fabsA[i] = (int)(fabsA[i] * delta + 0.45); //����
		if (fabsA[i])
			j++;
		else
			sgn[i] = 0;
	}
	if (!fabsA[len - 1])
		j++;
	temp = (int*)realloc(qtTHD.runs, j*sizeof(int));
	if (temp == NULL){
		printf("���·����ڴ�ʧ�ܣ�\n");
		exit(1);
	}
	qtTHD.runs = temp;
	qtTHD.lg = j;
	return qtTHD;
}

int quanTHD2(float *absA, float delta)  //ֻ��������
{
	extern int len;
	int i;
	for (i = 0; i < len; i++)
	{
		absA[i] = (int)(absA [i] * delta + 0.48); //����
	}
	return 0;
}