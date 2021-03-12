#ifndef _MODIFY_HIBERT2D_H_
#define _MODIFY_HIBERT2D_H_
#include<assert.h>

/*find_block.h*/
typedef struct
{
   int block[200];
   int order[200];
   int num;
}block;

typedef struct
{
   block basic_block;
   block h_order;
}total;


int Sgn(int value);										//判断input的正负
int power_of_2(int num);								//2的幂运算
block find_basic_block(int length);						//按2的幂次方分割成从大到小的数组
total judge_order(total basic, int w, int h);			//用于判断长度和宽度的计算顺序
block find_2s_inc(block bb, int dH);					//用于缩短要处理的数组长度
total judge_order_2(total basic, int w, int h);    //
block find_2s_inc_2(int H0, int dH, block inc0);



/*hilbert_wiki.h*/
int xy_point2d(int n, int x, int y);						//求出在特定大小图像中该xy坐标中对应的hilbert数
void rot(int n, int *x, int *y, int rx, int ry);
void d2xy(int n, int d, int *x, int *y);					//求出hilbert数在特定大小的图像中对应的xy坐标




/*draw_hilbert.h*/
int **create_arr(int width, int height);
int **draw_hibert(int W, int H);
int **HorizontalMirror(int **arr, int n, int m);
int *row(int n);
int **rotate_arr(int w, int h, int **old);




/*add.h*/
int **add(int **orginal_arr, int width, int height, int num);    //add
int **sub(int **orginal_arr, int width, int height, int num);    //sub

int *modify_Hilbert2D(int W, int H);


#endif
