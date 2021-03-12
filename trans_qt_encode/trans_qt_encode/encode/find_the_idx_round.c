#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

#define uint unsigned int
#define uchar unsigned char

/**
  * @brief
  * @param  n[IN] : 从索引0开始使用
  * @param  p[IN]
  * @param  thd_arr[IN]
  * @param  idx[OUT]
  * @retval none
  * @note
  */
void find_the_idx_round(Uint32_Dat* n, float p, int * thd_arr, int *idx)
{
	unsigned int sm = 0;
	float        psm = 0;
	int          lenn = 0;
	unsigned int ssm = 0;
	unsigned int sum = 0, sum2 = 0;

	/* sm=sum(n); */
	for (int i = 0; i < n->len; i++) {
		sm += n->dat[i];
		lenn++;
	}
	/* psm=p*sm; */
	psm = p * sm;
	/* idx=1; */
	*idx = 1;
	/* ssm=n(1); */
	ssm = n->dat[0];
	/* lenn=length(n); */
	lenn = lenn;

	/* 调试代码段 */
#if 0
	printf("psm  = %f\n", psm);
	printf("ssm  = %d\n", ssm);
	printf("lenn = %d\n", lenn);
#endif
	/* Hint:
	 C语言和MATLAB并没有实现索引同步，C语言指针索引减1 */
	//while (ssm < psm && ((*idx) - 1) < 24 && thd_arr[(*idx)] < lenn) {
	while (ssm < psm && ((*idx)) < 23 && thd_arr[(*idx)] < lenn) {
		/* sum(n(thd_arr(idx)+1:thd_arr(idx+1))); */
		sum = 0;
		for (int i = thd_arr[(*idx) - 1]; i < thd_arr[(*idx)]; i++) {
			sum += n->dat[i];
		}
		ssm += sum;
		//printf("ssm = %d sum = %d\n", ssm, sum);
		(*idx)++;
	}

	/* 调试代码段 */
#if 0
	printf("ssm=%d\n", ssm);
	printf("idx=%d\n", (*idx));
	while (1);
#endif

	if ((*idx) > 1) {
		/* sum(n(1:thd_arr(idx-1))) */
		sum = 0;
		for (int i = 0; i < thd_arr[(*idx) - 2]; i++)
			sum += n->dat[i];
		/* sum( n(1:thd_arr(idx))) */
		sum2 = 0;
		for (int i = 0; i < thd_arr[(*idx) - 1]; i++)
			sum2 += n->dat[i];
		/* if psm-sum(n(1:thd_arr(idx-1)))  <  sum( n(1:thd_arr(idx)))-psm */
		//if (psm - sum < sum2 - psm) {
		if ((psm - sum < sum2 - psm) || (lenn - (thd_arr[(*idx) - 1])) <= 1) 
		{
			(*idx)--;
		}
	}
	//STDOUT_INFO("idx=%d\n", *idx);

	return;
}

