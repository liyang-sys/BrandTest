//#include <stdio.h>
//#include "all.h"
//#include "parameter_setting.h"
//#include "encoding.h"
//
//float probility_estimation(Int32_Dat*rk, Int32_Dat*nk)
//{
//	float prob = 0;
//	float prb2 = (float)(rk->len) / (float)(sum(rk->dat, rk->len));
//	int thd = rk->dat[0];
//	for (int i = 0; i < rk->len; i++)
//	{
//		if (thd < (rk->dat[i]))
//		{
//			thd = (rk->dat[i]);
//		}
//	}
//	if (thd <= 7 && (rk->len) > 360)
//	{
//		prob = prob_est(nk);
//	}
//	else if(thd <= 12)
//	{
//		float prb1 = prob_est(nk);
//		prob = (prb1 + prb2) / 2.0;
//	}
//	else
//	{
//		prob = prb2;
//	}
//	return prob;
//}
#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char

void probility_estimation(Uint32_Dat* rk, Uint32_Dat* nk, float *prob)
{
	int sum_rk = 0;
	int max_rk = 0;
	float prb1 = 0;
	float prb2 = 0;
	int thd = 0;


	sum_rk = 0;
	for (int i = 0; i < rk->len; i++) {
		if (max_rk < rk->dat[i])  max_rk = rk->dat[i];
		sum_rk += rk->dat[i];
		//printf("%d:%d sum_rk = %d \n", i+1, rk->dat[i], sum_rk);
	}
	/* prb2=length(rk)/sum(rk); */
	prb2 = (float)rk->len / sum_rk;
	//printf("prb2 = %f\n", prb2);
	/* thd=max(rk); */
	thd = max_rk;
	/* if thd <=7 && length(rk)>360 */
	if (thd <= 7 && rk->len > 360) {
		/* prb1=prob_est(nk); */
		prob_est(nk, &prb1);
		/* prob=prb1; */
		*prob = prb1;
	}
	else if (thd <= 12) {  /* elseif thd<=12 */
	   /* prb1=prob_est(nk);  */
		prob_est(nk, &prb1);
		/* prob=(prb1+prb2)/2; */
		*prob = (prb1 + prb2) / 2;
	}
	else {
		/* prob=prb2; */
		*prob = prb2;
	}

	return;
}
