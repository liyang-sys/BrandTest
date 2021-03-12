//#include <stdio.h>
//#include "all.h"
//#include "encoding.h"
//
//#define uint unsigned int
//#define uchar unsigned char
//
//
///**
// * @brief
// * @param  z[IN]         : 空类型，根据sep_flg决定其类型，为Uint32_Dat或者Uint8_dat,会被销毁，请提前保存
// * @param  cr[OUT]
// * @param  cr_count[OUT] : 若不需要此变量，传递NULL
// * @param  r[OUT]        : 若不需要此变量，传递NULL
// * @param  _1B[OUT]      : 若不需要此变量，传递NULL
// * @param  p[OUT]        : 若不需要此变量，传递NULL
// * @param  nc[OUT]       : 若不需要此变量，传递NULL
// * @param  k[OUT]        : 若不需要此变量，传递NULL //李旸注释：输出比实际长度多1个，数组第一个的时候没有数据
// * @param  sep_flg[IN]   : 传入1代表为sep 0 1二进制码流，传入0代表是r数组
// * @retval none
// * @note
// */
////void K_Criterion(Uint8_Dat* z, float* cr, int* cr_count, Uint32_Dat* r, int* lB, float* p, Uint32_Dat* nc, int* k, int flg)
//void k_criterion(void* z, float* cr, int* cr_count, Uint32_Dat* r, int* lB, float* p, Uint32_Dat* nc, int* k, int sep_flg)
//{
//	float q = 0;
//	int dn = 0;
//	float an = 0;
//	int N0 = 0;
//	GOL gol;
//	int sum = 0;
//	int max = 0;
//
//	//int         k_temp = 0;
//	//Uint32_Dat  nc_temp;
//	//float       p_temp = 0;
//	//Uint32_Dat  r_temp;
//	//int         cr_count_temp = 0;
//	///* 如果不需要这些参数的返回值，则传参传NULL */
//	//if (k == NULL)          k = &k_temp;
//	//if (nc == NULL)         nc = &nc_temp;
//	//if (p == NULL)          p = &p_temp;
//	//if (r == NULL)          r = &r_temp;
//	//if (cr_count == NULL)   cr_count = &cr_count_temp;
//	int         k_temp = 0;
//	Uint32_Dat  nc_temp;
//	float       p_temp = 0;
//	Uint32_Dat  r_temp;
//	int         cr_count_temp = 0;
//	int         lB_temp = 0;
//	/* 如果不需要这些参数的返回值，则传参传NULL */
//	if (k == NULL)          k = &k_temp;
//	if (nc == NULL)         nc = &nc_temp;
//	if (p == NULL)          p = &p_temp;
//	if (r == NULL)          r = &r_temp;
//	if (cr_count == NULL)   cr_count = &cr_count_temp;
//	if (lB == NULL)         lB = &lB_temp;
//
//	/* if max(z)<=1 */
//	if (sep_flg == 1) {
//		/* [r,lB,p]=Golomb0(z) ; */
//		gol = Golomb0(((Uint8_Dat *)z)->dat, ((Uint8_Dat *)z)->len);
//		r->dat = gol.r;
//		r->len = gol.lenr;
//		*lB = gol.lastBit;
//		*p = gol.p;
//		/* N0=length(r); */
//		N0 = r->len;
//		//修改的
//		//free(((Uint8_Dat *)z)->dat);
//		//((Uint8_Dat *)z)->dat = NULL;
//	}
//	else {
//		/* r=z; */
//		r->dat = ((Uint32_Dat *)z)->dat;
//		r->len = ((Uint32_Dat *)z)->len;
//		//((Uint32_Dat *)z)->dat = NULL; 修改的
//		/* N0=length(r); */
//		N0 = r->len;
//		/* sum(r); */
//		sum = 0;
//		for (int i = 0; i < r->len; i++) sum += r->dat[i];
//		*p = (float)N0 / sum;
//		*lB = 0;
//	}
//
//	/* q=1-p; */
//	q = 1 - *p;
//
//	/* 一些辅助变量，因为log函数需要的时间比较久，所以有些固定参数的值就求出来存放在一遍就行了 */
//	/* 减少调用log函数的次数 */
//	float log_q = log(q);
//
//	/* max(r) */
//	max = 0;
//	for (int i = 0; i < r->len; i++) {
//		if (max < r->dat[i]) max = r->dat[i];
//	}
//
//	/* nc=hist(r, [1:max(r)]); */
//	nc->len = max + 1;
//	nc->dat = (unsigned int *)malloc(sizeof(unsigned int)*(nc->len + 1));
//	memset(nc->dat, 0, sizeof(unsigned int)*(nc->len + 1));
//	for (int i = 0; i < r->len; i++) {
//		nc->dat[(unsigned int)(r->dat[i])] += 1;
//	}
//
//	/* 调试代码段 */
//#if 0
//	int sum = 0;
//	STDOUT_DEBUG("nc数组输出\n");
//	for (int i = 0; i < nc->len; i++) {
//		STDOUT_DEBUG("%4d:%2d\n", i + 1, nc->dat[i]);
//		sum += nc->dat[i];
//	}
//	STDOUT_DEBUG("长度：%d 和为%d\n", nc->len, sum);
//	while (1);
//#endif
//
//	/* k=floor(-log(N0)/log(q)); */
//	*k = floor(-log(N0) / log_q);
//
//	/* cr_count=sum(nc(k+1:length(nc))); */
//	*cr_count = 0;
//	for (int i = (*k) + 1; i < nc->len; i++)
//		*cr_count += nc->dat[i];
//
//	/* cr=0; dn=0; */
//	*cr = 0; dn = 0;
//	for (int n = (*k) + 1; n < nc->len; n++) {
//		/* if nc(n)>0 */
//		if (nc->dat[n] > 0) {
//			/* if dn==0 */
//			if (dn == 0) {
//				/* an=nc(n)/N0; */
//				an = (float)nc->dat[n] / N0;
//				/* cr = cr + log(an/p) -(n-1)*log(q) ; */
//				*cr += log(an / (*p)) - (n - 1)*log_q;
//				//STDOUT_DEBUG("(1)cr = %f\n", *cr);
//			}
//			else {
//				/* an=(nc(n)/N0)/(dn+1);   */
//				an = ((float)nc->dat[n] / N0) / (dn + 1);
//				/*  cr = cr + log(an/p) -(n-1-ceil(dn/2))*log(q) ; */
//				*cr += log(an / (*p)) - (n - 1 - ceil((float)dn / 2))*log_q;
//				//STDOUT_DEBUG("(2)cr = %f dn=%d log(q)=%f\n", *cr, dn, log_q);
//			}
//			/* dn=0; */
//			dn = 0;
//		}
//		else {
//			/* dn=0; */
//			dn += 1;
//		}
//	}
//
//	/* 调试代码段 */
//#if 0
//	STDOUT_DEBUG("cr = %f\n", *cr);
//	STDOUT_DEBUG("cr_count = %d\n", *cr_count);
//#endif
//	if (nc == &nc_temp) {
//		free(nc->dat);
//	}
//	if (sep_flg == 1)
//	{
//		if (r == &r_temp) {
//			free(r->dat);
//		}
//	}
//
//	//修改的
//	//if (sep_flg == 1)
//	//{
//	//	free(r->dat);
//	//}
//	return;
//}

#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char


/**
 * @brief
 * @param  z[IN]         : 空类型，根据sep_flg决定其类型，为Uint32_Dat或者Uint8_dat,会被销毁，请提前保存
 * @param  cr[OUT]       : 若不需要此变量，传递NULL
 * @param  cr_count[OUT] : 若不需要此变量，传递NULL
 * @param  r[OUT]        : 若不需要此变量，传递NULL
 * @param  _1B[OUT]      : 若不需要此变量，传递NULL
 * @param  p[OUT]        : 若不需要此变量，传递NULL
 * @param  nc[OUT]       : 若不需要此变量，传递NULL, 实际使用索引从1开始
 * @param  k[OUT]        : 若不需要此变量，传递NULL //李旸注释：输出长度没变数据多一个，数组第一个的时候没有数据
 * @param  sep_flg[IN]   : 传入1代表为sep 0 1二进制码流，传入0代表是r数组
 * @retval none
 * @note
 */
void k_criterion(void* z, float* cr, int* cr_count, Uint32_Dat* r, int* lB, float* p, Uint32_Dat* nc, int* k, int sep_flg)
{
	float   q = 0;
	int     dn = 0;
	float   an = 0;
	int     N0 = 0;
	GOL     gol;
	int     sum = 0;
	int     max = 0;

	int         k_temp = 0;
	Uint32_Dat  nc_temp;
	float       p_temp = 0;
	Uint32_Dat  r_temp;
	int         cr_count_temp = 0;
	int         lB_temp = 0;
	float       cr_temp;
	/* 如果不需要这些参数的返回值，则传参传NULL */
	if (cr == NULL)         cr = &cr_temp;
	if (k == NULL)          k = &k_temp;
	if (nc == NULL)         nc = &nc_temp;
	if (p == NULL)          p = &p_temp;
	if (r == NULL)          r = &r_temp;
	if (cr_count == NULL)   cr_count = &cr_count_temp;
	if (lB == NULL)         lB = &lB_temp;

	/* if max(z)<=1 */
	if (sep_flg == 1) {
		/* [r,lB,p]=Golomb0(z) ; */
		gol = Golomb0(((Uint8_Dat *)z)->dat, ((Uint8_Dat *)z)->len);
		r->dat = gol.r;
		r->len = gol.lenr;
		*lB = gol.lastBit;
		*p = gol.p;
		/* N0=length(r); */
		N0 = r->len;
	}
	else {
		/* r=z; */
		r->dat = ((Uint32_Dat *)z)->dat;
		r->len = ((Uint32_Dat *)z)->len;
		/* N0=length(r); */
		N0 = r->len;
		/* sum(r); */
		sum = 0;
		for (int i = 0; i < r->len; i++) sum += r->dat[i];
		*p = (float)N0 / sum;
		*lB = 0;
	}

	/* q=1-p; */
	q = 1 - *p;

	/* 一些辅助变量，因为log函数需要的时间比较久，所以有些固定参数的值就求出来存放在一遍就行了 */
	/* 减少调用log函数的次数 */
	float log_q = log(q);

	/* max(r) */
	max = 0;
	for (int i = 0; i < r->len; i++) {
		if (max < r->dat[i]) max = r->dat[i];
	}

	/* nc=hist(r, [1:max(r)]); */
	nc->len = max;
	nc->dat = (unsigned int *)malloc(sizeof(unsigned int)*(nc->len + 1));
	memset(nc->dat, 0, sizeof(unsigned int)*(nc->len + 1));
	for (int i = 0; i < r->len; i++) {
		nc->dat[(unsigned int)(r->dat[i])] += 1;
	}

	/* 调试代码段 */
#if 0
	int sum = 0;
	STDOUT_DEBUG("nc数组输出\n");
	for (int i = 1; i <= nc->len; i++) {
		STDOUT_DEBUG("%4d:%2d\n", i + 1, nc->dat[i]);
		sum += nc->dat[i];
	}
	STDOUT_DEBUG("长度：%d 和为%d\n", nc->len, sum);
	while (1);
#endif

	/* k=floor(-log(N0)/log(q)); */
	*k = floor(-log(N0) / log_q);

	/* cr_count=sum(nc(k+1:length(nc))); */
	*cr_count = 0;
	for (int i = (*k) + 1; i <= nc->len; i++)
		*cr_count += nc->dat[i];

	/* cr=0; dn=0; */
	*cr = 0; dn = 0;
	for (int n = (*k) + 1; n <= nc->len; n++) {
		/* if nc(n)>0 */
		if (nc->dat[n] > 0) {
			/* if dn==0 */
			if (dn == 0) {
				/* an=nc(n)/N0; */
				an = (float)nc->dat[n] / N0;
				/* cr = cr + log(an/p) -(n-1)*log(q) ; */
				*cr += log(an / (*p)) - (n - 1)*log_q;
				//STDOUT_DEBUG("(1)cr = %f\n", *cr);
			}
			else {
				/* an=(nc(n)/N0)/(dn+1);   */
				an = ((float)nc->dat[n] / N0) / (dn + 1);
				/*  cr = cr + log(an/p) -(n-1-ceil(dn/2))*log(q) ; */
				*cr += log(an / (*p)) - (n - 1 - ceil((float)dn / 2))*log_q;
				//STDOUT_DEBUG("(2)cr = %f dn=%d log(q)=%f\n", *cr, dn, log_q);
			}
			/* dn=0; */
			dn = 0;
		}
		else {
			/* dn=0; */
			dn += 1;
		}
	}

	/* 调试代码段 */
#if 0
	STDOUT_DEBUG("cr = %f\n", *cr);
	STDOUT_DEBUG("cr_count = %d\n", *cr_count);
#endif
	if (nc == &nc_temp) {
		free(nc->dat);
	}


	return;
}

