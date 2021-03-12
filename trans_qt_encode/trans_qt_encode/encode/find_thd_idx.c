#include <stdio.h>
#include <math.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

int  find_thd_idx(Uint32_Dat*n,float p,int *th) //n的真实数据得从1开始，序号等于零时置零
{
	int sm = 0;
	for (int i = 0; i < n->len; i++)
	{
		sm += n->dat[i];
	}
	float psm = p * sm;
	int idx = 1;
	int ssm = n->dat[0];
	while (ssm < psm)
	{
		for (int i = th[idx] + 1; i <= th[idx + 1]; i++)
		{
			ssm += n->dat[i-1];
		}
		idx++;
	}
	return idx;
}