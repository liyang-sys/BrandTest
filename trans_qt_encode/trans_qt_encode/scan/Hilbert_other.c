#include "modify_Hilbert3D.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

static long int p_2s[14] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,8192 };
//计算一个数可以拆分为一个用2的N次方组成的数组，结果均为正数
block0 find_basic_block0(int length)
{
	block0 b_b;
	int i, j;
	memset(&b_b, 0, sizeof(b_b));    //将block b_b置零
	for (i = 0; length > 0; i++){
		for (j = 0; p_2s[j] <= length; j++){
			b_b.block[i] = p_2s[j];
		}
		b_b.order[i] = j - 1;
		length = length - b_b.block[i];
	}
	b_b.num = i;
	return (b_b);
}

//用于判断符号
int Sign(int value)
{
	if (value > 0)
		return 1;
	else if(value < 0)
		return -1;
	else 
		return 0;
}


block0 find_2s_inc0(int H0, int dH, block0 inc0)    //缩短分解后的序列
{
	int n_inc = inc0.num;
	int *inc = row0(inc0.num);
	int ia = 1, i, j;
	block0 dh;
	for (i = 0; i < inc0.num; i++)
		inc[i] = inc0.block[i];
	while (ia < inc0.num - 1)
	{
		if (inc0.block[ia - 1] / 2 == inc0.block[ia] && abs(dH)>abs(inc0.block[ia - 1] + inc0.block[ia]))
		{
			inc[ia - 1] = inc0.block[ia - 1] * 2;
			dH = dH - inc[ia - 1];
			dh = find_basic_block0(abs(dH));
			n_inc = ia + dh.num;
			ia++;
			for (i = ia - 1, j = 0; i < n_inc; j++, i++)
				inc[i] = Sign(dH)*dh.block[j];
			for (i = 0; i < inc0.num; i++)
			if (i < n_inc)
				inc0.block[i] = inc[i];
			else
				inc0.block[i] = 0;
			inc0.num = n_inc;
		}
		else
		{
			dH = dH - inc[ia - 1];
			ia++;
		}
	}
	return inc0;
}


/*用于动态创建任意长度二维数组并返回*/
int **create_arr0(int width, int height)
{
	int **new_arr = (int **)malloc(width*sizeof(int *));										//创建行
   for (int i = 0; i < width; i++)
   {
      new_arr[i] = (int *)malloc(height*sizeof(int));											//创建列
      assert(new_arr[i] != NULL);
   }
	return new_arr;
}

int *row0(int n)
{
	int *arr = (int *)malloc(n*sizeof(int *));
   assert(arr != NULL);
	return arr;
}

//改变finish的值
void swap(int *a, int *b)
{
	int temp = *a;
	*a = *b;
	*b = temp;
}

//三维的按行扫描
int **zigzag3D(int Nr, int Nc, int rank)
 {
	int i, j, ii;
	int Nd = pow(2, rank);
	int *dp = row0(Nc), *dn = row0(Nc);
	int **m = create_arr0(3, Nr*Nc*Nd);
	int pos = 0, dir = 1, dirfm = 1;
	int d;
	int *row1 = row0(Nr);
	int r;

	for (i = 0; i < Nc; i++)
		dp[i] = i+1;
	for (i = 0, j = Nc; i<Nc, j>0; i++, j--)
		dn[i] = j;
	for (i = 0; i < 3;i++)
	for (j = 0; j < Nr*Nc*Nd; j++)
		m[i][j] = 0;
	for (d = 1; d <= Nd; d++)
	{
		if (dirfm>0)
		for (i = 0; i < Nr; i++)
			row1[i] = i + 1;
		else
		for (i = 0, j = Nr; i<Nr, j>0; i++, j--)
			row1[i] = j;
		for (i = 0; i < Nr; i++)
		{
			ii = 0;
			r = row1[i];
			if (dir>0)
			for (j = pos; j < pos + Nc ; j++)
			{
				m[0][j] = r;
				m[2][j] = d;
				m[1][j] = dp[ii];
				ii++;
			}
			else
			for (j = pos; j < pos + Nc ; j++)
			{
				m[0][j] = r;
				m[2][j] = d;
				m[1][j] = dn[ii];
				ii++;
			}
			pos = pos + Nc;
			dir = -dir;
		}
		dirfm = -dirfm;
	}
	//释放空间													
	free(dp); free(dn); free(row1);

	return m;
}