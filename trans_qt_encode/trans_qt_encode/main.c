#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "quant_sub.h"
#include "transform.h"
#include <math.h>
#include "encoding.h"

unsigned char **sn;
unsigned char *snbin;
unsigned char *bin;
int ptr = 0;
int **nc;
int len;
float ***VideoData;
float ***reVideoData;
float ***PTVData;
float ***rePTVData;
float *absA2;
float *absA;
float *T;
float *T2;
union Fabs f1;

void memoryAllocation(int H, int W, int m, int n)
{
	if ((VideoData = (float***)malloc(sizeof(float**)* 96)) == NULL){
		printf("����VideoData�����һάʧ�ܣ�\n");
		exit(1);
	}
	for (int i = 0; i < 96; i++){
		if ((VideoData[i] = (float**)malloc(sizeof(float*)* H)) == NULL){
			printf("����VideoData����ڶ�άʧ�ܣ�\n");
			exit(1);
		}
		for (int j = 0; j < H; j++){
			if ((VideoData[i][j] = (float*)malloc(sizeof(float)* W)) == NULL){
				printf("����VideoData�������άʧ�ܣ�\n");
				exit(1);
			}
		}
	}
	if ((reVideoData = (float***)malloc(sizeof(float**)* 96)) == NULL){
		printf("����VideoData�����һάʧ�ܣ�\n");
		exit(1);
	}
	for (int i = 0; i < 96; i++){
		if ((reVideoData[i] = (float**)malloc(sizeof(float*)* H)) == NULL){
			printf("����VideoData����ڶ�άʧ�ܣ�\n");
			exit(1);
		}
		for (int j = 0; j < H; j++){
			if ((reVideoData[i][j] = (float*)malloc(sizeof(float)* W)) == NULL){
				printf("����VideoData�������άʧ�ܣ�\n");
				exit(1);
			}
		}
	}
	if ((PTVData = (float***)malloc(sizeof(float**)* 96)) == NULL){
		printf("����PTVData�����һάʧ�ܣ�\n");
		exit(1);
	}
	for (int i = 0; i < 96; i++){
		if ((PTVData[i] = (float**)malloc(sizeof(float*)* m)) == NULL){
			printf("����PTVData����ڶ�άʧ�ܣ�\n");
			exit(1);
		}
		for (int j = 0; j < m; j++){
			if ((PTVData[i][j] = (float*)malloc(sizeof(float)* n)) == NULL){
				printf("����PTVData�������άʧ�ܣ�\n");
				exit(1);
			}
		}
	}
	if ((rePTVData = (float***)malloc(sizeof(float**)* 96)) == NULL){
		printf("����PTVData�����һάʧ�ܣ�\n");
		exit(1);
	}
	for (int i = 0; i < 96; i++){
		if ((rePTVData[i] = (float**)malloc(sizeof(float*)* m)) == NULL){
			printf("����PTVData����ڶ�άʧ�ܣ�\n");
			exit(1);
		}
		for (int j = 0; j < m; j++){
			if ((rePTVData[i][j] = (float*)malloc(sizeof(float)* n)) == NULL){
				printf("����PTVData�������άʧ�ܣ�\n");
				exit(1);
			}
		}
	}
	if ((absA2 = (float*)malloc(sizeof(float)*len)) == NULL){
		printf("����absA2����ʧ�ܣ�\n");
		exit(1);
	}
	if ((absA = (float*)malloc(sizeof(float)*len)) == NULL){
		printf("����absA2����ʧ�ܣ�\n");
		exit(1);
	}
	if ((T = (float*)malloc(sizeof(float)* 131072)) == NULL){
		printf("����T����ʧ�ܣ�\n");
		exit(1);
	}
	if ((T2 = (float*)malloc(sizeof(float)* 131073)) == NULL){
		printf("����T2����ʧ�ܣ�\n");
		exit(1);
	}
	if ((f1.qcf = (float**)malloc(sizeof(float*)* 8)) == NULL){
		printf("����qcf�����һάʧ�ܣ�\n");
		exit(1);
	}
	for (int i = 0; i < 8; i++){
		if ((f1.qcf[i] = (float*)malloc(sizeof(float)*len)) == NULL){
			printf("����qcf����ڶ�άʧ�ܣ�\n");
			exit(1);
		}
	}
	if ((sn = (unsigned char**)calloc(8, sizeof(char*))) == NULL){
		printf("����sn����ʧ�ܣ�\n");
		exit(1);
	}
	for (int i = 0; i < 8; i++){
		if ((sn[i] = (unsigned char*)calloc(len, sizeof(char))) == NULL){
			printf("����sn����ڶ�άʧ�ܣ�\n");
			exit(1);
		}
	}
	if ((snbin = (unsigned char*)calloc((len >> 3) + 4, sizeof(char))) == NULL){
		printf("����snbin����ʧ�ܣ�\n");
		exit(1);
	}
	if ((bin = (unsigned char*)calloc(len + 4, sizeof(char))) == NULL){
		printf("����bin����ʧ�ܣ�\n");
		exit(1);
	}
	if ((nc = (int**)calloc(8, sizeof(int*))) == NULL){
		printf("����nc�����һάʧ�ܣ�\n");
		exit(1);
	}
	for (int i = 0; i < 8; i++){
		if ((nc[i] = (int*)calloc(len * 10, sizeof(int))) == NULL){
			printf("����nc����ڶ�άʧ�ܣ�\n");
			exit(1);
		}
	}
}

//int main()
//{
//	test();
//}
int main()
{
	int W, H, n, m;
	int *m0[38], w[7], h[7];
	float T0, delta;
	FILE *fp3;
	errno_t err;
	printf("Please input the width\n");   //W��
	scanf_s("%d", &W);
	while (W % 8 != 0)
	{
		printf("Your input is illegal,Please try again\n");
		scanf_s("%d", &W);
	}

	printf("Please input the height\n");   //H��
	scanf_s("%d", &H);
	while (H % 8 != 0)
	{
		printf("Your input is illegal,Please try again\n");
		scanf_s("%d", &H);
	}

	printf("������delta��\n");
	scanf_s("%f", &delta, sizeof(float));
	m = H / 8;
	n = W / 8;
	len = n*m * 8;
	printf("׼��������ʼ:\n");
	//�ڴ����
	memoryAllocation(H, W, m, n);
	//��ȡ�������
	printf("��ʼ���ݶ�ȡ��\n");

	if ((err = fopen_s(&fp3, "Tcoef.dat", "rb")) != 0){
		printf("��T����ʧ�ܣ�\n");
		exit(1);
	}
	fseek(fp3, 0, SEEK_SET);
	fread(T, sizeof(float), 131071, fp3);
	fclose(fp3);
	printf("���ݶ�ȡ��ɣ�\n");

	T0 = 0.55*delta;
	T2[1] = T0*T0; T2[0] = 0; //�����T��T2�ĳ��ȶ�������һ��������Ϊ��matlabƥ�䣬��1��ʼ������T[0]�൱��Ϊ��
	for (int i = 131070; i >= 0; i--){
		T[i + 1] = T[i] * delta;
		T2[i + 2] = T[i+1] * T[i+1];
	}
	printf("׼��������ɣ�\n\n");
	//.....���任��ʼ��.....
	clock_t st2, end2;
	unsigned char*** imageData = forward_transform(W, H); //VideoData��PTV����δ�����

	//���Բ���
	st2 = clock();
	GetM(W, H, w, h, m0);
	printf("\nHilbert·����ȡ��ϣ�\n\n");
	printf(".....������ʼ��.....\n");
	

	/*��һ��32*/
	/*��Ҫ��harr�Ĳ���*/
	quant_subDC_harr(0, w, h, m0, delta);//DC
	quant_sub7_harr(0, w, h, m0, delta);//DC���7���Ӵ�
	for (int i = 1; i < 8; i++)
		quant_sub8_harr(0, i, w, h, m0, delta);//8���Ӵ�
	///*����Ҫharr����*/
	for (int j = 1; j < 4; j++){
		quant_subDC_noharr(j, w, h, m0, delta);//����Ҫ��harr��DC
		quant_sub7_noharr(j, w, h, m0, delta);//DC���7���Ӵ�
		for (int i = 1; i < 8; i++)
			quant_sub8_noharr(j, i, w, h, m0, delta);//8���Ӵ�
	}
	/*�ڶ���32*/
	quant_subDC_harr(32, w, h, m0, delta);//DC
	quant_sub7_harr(32, w, h, m0, delta);//DC���7���Ӵ�
	for (int i = 1; i < 8; i++)
		quant_sub8_harr(32, i, w, h, m0, delta);//8���Ӵ�
	///*����Ҫharr����*/
	for (int j = 33; j < 36; j++){
		quant_subDC_noharr(j, w, h, m0, delta);//����Ҫ��harr��DC
		quant_sub7_noharr(j, w, h, m0, delta);//DC���7���Ӵ�
		for (int i = 1; i < 8; i++)
			quant_sub8_noharr(j, i, w, h, m0, delta);//8���Ӵ�
	}
	/*������32*/
	quant_subDC_harr(64, w, h, m0, delta);//DC
	quant_sub7_harr(64, w, h, m0, delta);//DC���7���Ӵ�
	for (int i = 1; i < 8; i++)
		quant_sub8_harr(64, i, w, h, m0, delta);//8���Ӵ�
	///*����Ҫharr����*/
	for (int j = 65; j < 68; j++){
		quant_subDC_noharr(j, w, h, m0, delta);//����Ҫ��harr��DC
		quant_sub7_noharr(j, w, h, m0, delta);//DC���7���Ӵ�
		for (int i = 1; i < 8; i++)
			quant_sub8_noharr(j, i, w, h, m0, delta);//8���Ӵ�
	}
	end2 = clock();
	double duration2 = (double)(end2 - st2) / CLOCKS_PER_SEC;
	printf("\n\nʱ�� %f seconds\n", duration2);
	printf("%f\n", rePTVData[3][0][0]);
	for (int i = 0; i < 8; i++){
		for (int j = 0; j < 8; j++)
			printf("%f ", reVideoData[3][i][j]);
		printf("\n");
	}
	printf("\n\n");
	//rearrange
#if 0
	FILE *fp2;
	if ((err = fopen_s(&fp2, "E:\\�������\\����\\2021-01-15\\test_en_sub3d_sub2_sub\\PTVoutput_rearrange.dat", "wb")) == 0){
		for (int i = 0; i < 96; i++)
		{
			for (int j = 0; j < m; j++)
			{
				fwrite(&PTVData3[i][j][0], sizeof(float), n, fp2);
			}
		}
		fclose(fp2);
	}
#endif
	/*--------------------------------------------------------------------------------*/
#if 1
	//.....��任��ʼ��.....
	inverse_transformation(W, H);

	for (int i = 0; i < 8; i++){
		for (int j = 0; j < 8; j++)
			printf("%f ", reVideoData[0][i][j]);
		printf("\n");
	}
	printf("\n\n");
	//if ((err = fopen_s(&fp3, "E:\\�������\\����\\2021-01-15\\test_en_sub3d_sub2_sub\\VideoData.dat", "wb")) != 0){
	//	printf("������ʧ�ܣ�\n");
	//	system("pause");
	//	return 0;
	//}
	//for (int i = 0; i < 96; i++)
	//{
	//	for (int j = 0; j<1080; j++)
	//	{
	//		fwrite(&reVideoData[i][j][0], sizeof(float), 1920, fp3);
	//	}
	//}
	//fclose(fp3);
#endif
	/*-----------------------------------�������---------------------------------------------*/
	float snr;
	float Npix = (float)H*W;
	float e, sume = 0.0;
	float sumsnr = 0, minsnr = 100.0;
	for (int i = 0; i < 96; i++){
		sume = 0.0;
		for (int j = 0; j < 1080; j++){
			for (int k = 0; k < 1920; k++){
				e = (float)imageData[i][j][k] - reVideoData[i][j][k];
				e *= e;
				sume += e;
			}
		}
		float er = sqrt((double)(sume / Npix));
		snr = 20.0*log10((double)(255.0 / er));
		//printf("%f\n", snr);
		sumsnr += snr;
		if (snr < minsnr)
			minsnr = snr;
	}
	float meansnr = sumsnr / 96.0;

	printf("mean(snr)=%f , min(snr)=%f\n", meansnr, minsnr);

	printf("���������\n\n");
	system("pause");
	return 0;
}

