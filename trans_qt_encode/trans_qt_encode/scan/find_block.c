#include "modify_Hilbert2D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

long int p_2s[14] = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192};

//用于判断正负
int Sgn(int value)
{
	if (value >= 0)
		return 1;
	else
		return -1;
}

//用于计算2的幂次方
int power_of_2(int num)
{
	int i, result = 1;
	if (num < 0)
		return 0;
	for (i = 0; i < num; i++)
		result = result * 2;
	return result;
}


//此函数与MATLAB中find_2spower功能一直，把input拆分成2的幂次方的组合
block find_basic_block(int length)
{
	block b_b;
	int i, j;
	memset(&b_b, 0, sizeof(b_b));
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

//用于判断将要计算的长度和宽度的顺序
total judge_order(total basic, int w, int h)
{
   int i, j;
   int H = h;
   int dH;
   int count;
   int sum = 0;
   int time = 0;
   int sgn;
   int power;
   int length;
   int arr_w[200] = { 0 }, arr_h[200] = { 0 };
   block bb_w, bb_h;
   block h_new;
   bb_w = find_basic_block(w);												//求出宽度的分割数组
   bb_h = find_basic_block(h);												//求出高度的分割数组
   length = bb_w.block[0];

   if (h >= bb_w.block[0]){												//若满足此要求，可直接传出bb_w和bb_h
      basic.basic_block = bb_w;
      basic.h_order = bb_h;
   }
   //如果h<bb.block[0],有两种情况，若h<3/4bb.block[0],使用加，若h>3/4bb.block[0]，使用减
   else if (h < bb_w.block[0]){
      for (i = 0; bb_w.block[i] > h; i++){
         length = bb_w.block[i];
         for (count = 0; h < length; count++)
            length = length >> 1;
         if (h >= length && h < length + (length << 1) / 4)
            power = power_of_2(count);
         else if (h < (length << 1) && h >= length + (length << 1) / 4)
            power = power_of_2(count - 1);
         sum += power;
         time++;
      }
      ///第一种
      if (h >= length && h < length + (length << 1) / 4){
         for (i = 0; bb_h.block[i] != 0; i++)
            arr_h[i] = bb_h.block[i];									//h就不用变，w序列也许会变

         for (i = sum, j = time; bb_w.block[j] != 0; i++, j++)			//全体数后移count*2的位置
            arr_w[i] = bb_w.block[j];
         for (i = 0; i < sum; i++)										//使位移后前面的数等于算出来的length，此时的length刚好<h
            arr_w[i] = length;
      }
      ////第二种
      else if (h < (length << 1) && h >= length + (length << 1) / 4){		//h和w都变了
         length = length << 1;											//此时length刚好>h
         dH = h - length;
         sgn = Sgn(dH);
         h_new = find_basic_block(abs(dH));
         arr_h[0] = length;
         for (i = 1; i <= h_new.num; i++)								//此处是改变高度
            arr_h[i] = sgn*h_new.block[i - 1];

         for (i = sum, j = time; bb_w.block[j] != 0; i++, j++)
            arr_w[i] = bb_w.block[j];

         for (i = 0; i < sum; i++)										//这里是改变宽度
            arr_w[i] = length;
      }
      for (i = 0; arr_w[i] != 0; i++)										//把计算好的序列送出去
         basic.basic_block.block[i] = arr_w[i];

      i = 0;
      while (basic.basic_block.block[i] != 0){
         for (int j = 0; j < 14; j++){
            if (basic.basic_block.block[i] == p_2s[j]){
               basic.basic_block.order[i] = j;
               i++;
            }
         }
      }
      basic.basic_block.num += sum - 1;
      for (i = 0; arr_h[i] != 0; i++)
         basic.h_order.block[i] = arr_h[i];
      for (; i <= basic.h_order.num; i++)
         basic.h_order.block[i] = 0;

      i = 0;
      int t;
      basic.h_order.num = 0;
      while (basic.h_order.block[i] != 0){
         if (basic.h_order.block[i] > 0) t = basic.h_order.block[i];
         else if (basic.h_order.block[i] < 0) t = -basic.h_order.block[i];
         for (int j = 0; j < 14; j++){
            if (t == p_2s[j]){
               basic.h_order.order[i] = j;
               i++;
               basic.h_order.num++;
            }
         }
      }

   }
   return basic;
}

//此函数与MATLAB中find_2sinc功能一致，用于判断能否缩短所需计算数组
block find_2s_inc(block bb,int dH)
{
	int i, j, ii;
	block b_b;
	int ia;
	int H = bb.block[0];
	int *inc;
	int len = bb.num;
	inc = (int *)malloc((len*2) * sizeof(int));								//动态建立一个用来计算新的组合的数组
	
	for (i = 0; i < len*2; i++)
		inc[i] = 0;
	for (i = 0; i < bb.num; i++)											//使该数组等于原有数组
		inc[i] = bb.block[i + 1];
	for (ia = 0; ia < bb.num - 1;){											//从0开始遍历到n-1个
		if (inc[ia] >> 1 == inc[ia + 1] && abs(dH) > abs(inc[ia] + inc[ia + 1])){    
			inc[ia] = inc[ia] << 1;
			dH = dH - inc[ia];
			b_b = find_basic_block(abs(dH));
			bb.num = ia + b_b.num;
			ia++;
			for (j = ia, ii = 0; j <= bb.num; ii++, j++)
				inc[j] = Sgn(dH)*b_b.block[ii];
			for (j = bb.num + 1; j < len + 1; j++)
				inc[j] = 0;
		}
		else{
			dH = dH - inc[ia];
			ia = ia + 1;
		}
	}
	for (i = 0; i <= bb.num; i++){
		bb.block[i + 1] = inc[i];
	}
	for (i += 1; i < len; i++)
		bb.block[i] = 0;

	free(inc);
	return bb;
}





//x是横坐标，y是纵坐标，分别对应图像的二维数组坐标，n是阶次
//此函数作为画出N阶次下的希尔伯特曲线
int xy_point2d(int n, int x, int y)
{
   int rx, ry, s, d = 0;
   for (s = n / 2; s>0; s /= 2) {
      //由于s是2的某次方，即在该位为1，而其他位数全为0，
      //则如果x<s，x&s=0，x>=s ,x&s>0
      rx = (x & s) > 0;
      ry = (y & s) > 0;
      d += s * s * ((3 * rx) ^ ry);
      rot(s, &x, &y, rx, ry);
   }
   return d;
}

//convert d to (x,y)
void d2xy(int n, int d, int *x, int *y)
{
   int rx, ry, s, t = d;
   *x = *y = 0;
   for (s = 1; s<n; s *= 2) {
      rx = 1 & (t / 2);
      ry = 1 & (t ^ rx);
      rot(s, x, y, rx, ry);
      *x += s * rx;
      *y += s * ry;
      t /= 4;
   }
}

//旋转象限
void rot(int n, int *x, int *y, int rx, int ry)
{
   if (ry == 0) {
      if (rx == 1) {
         *x = n - 1 - *x;
         *y = n - 1 - *y;
      }
      //Swap x and y
      int t = *x;
      *x = *y;
      *y = t;
   }
}

//新写的2021-1-15，对find_basic_block得出的对w,和h分解的系数做进一步处理。
total judge_order_2(total basic, int w, int h)
{
   block t;
   int delH;
   int tnum = basic.h_order.num;

   if (h >= basic.basic_block.block[0]){												//若满足此要求，可直接传出bb_w和bb_h
      for (int i = 1; i <= basic.h_order.num; i++){
         basic.h_order.block[i - 1] = basic.h_order.block[i];
      }
      basic.h_order.num--;
      return basic;
   }
   delH = h - basic.basic_block.block[0];
   if (delH != 0){
      if (delH > 0){
         basic.h_order = find_basic_block(delH);
      }
      else if (delH < 0){
         t = find_basic_block(-delH);
         for (int i = 0; i < t.num; i++){
            basic.h_order.block[i] = -t.block[i];
            basic.h_order.order[i] = t.order[i];
         }
         basic.h_order.num = t.num;
      }      
   }
   else{
      basic.h_order.num = 0;
   }

   for (int i = basic.h_order.num; i < tnum; i++){
      basic.h_order.block[i] = 0;
   }

   return basic;
}


//find_2s_inc的重新实现
//block find_2s_inc_2(block h, int dH)
//{
//   int H0 = h.block[0];
//   int n_inc = h.num;
//   int ia = 2;
//   int t;
//   int sum = 0;
//   block H;
//
//   if (dH > 0)t = dH;
//   else t = -dH;
//   for (int i = 1; i <= ia - 1; i++){
//      sum += h.block[i];
//   }
//   while (ia <= n_inc && (H0 <= (t << 1) || sum + h.block[ia] << 1 + H0 <= 0))
//   {
//      dH = dH - h.block[ia];
//      if (dH > 0)t = dH;
//      else t = -dH;
//      ia++;
//      sum += ia - 1;
//   }
//
//   int i = 0;
//   while (h.block[i] != 0){
//      H.block[i] = h.block[i];
//      i++;
//   }
//
//   if (dH > 0)t = dH;
//   else t = -dH;
//   int tHi;
//   if (H.block[ia] > 0)tHi = H.block[ia];
//   else tHi = -H.block[ia];
//   while  (ia<n_inc - 1)
//   if (H.block[ia] >> 1 == H.block[ia + 1] && t>tHi + H.block[ia + 1]){
//      h.block[ia] = bitshift(inc0(ia), 1);
//      dH = dH - inc(ia);
//      [idx, h, nh] = find_2s_power(abs(dH), p2s); %%%%%%%%%%%%%%%
//         n_inc = ia + nh; ia = ia + 1;
//      inc(ia:n_inc) = sign(dH)*h;  inc = inc(1:n_inc); inc0 = inc;
//   }
//
//   else{
//      dH = dH - inc(ia);
//      ia = ia + 1;
//   }
//
//
//}

block find_2s_inc_2(int H0, int dH, block inc0)    //缩短分解后的序列
{
   int n_inc = inc0.num;
   int *inc = (int *)malloc(n_inc*sizeof(int));
   int ia = 1, i, j;
   block dh;
   int t, sum = 0;
   for (i = 0; i < inc0.num; i++)
      inc[i] = inc0.block[i];

   if (dH > 0)t = dH;
      else t = -dH;
   for (int i = 0; i < ia - 1; i++){
      sum += inc0.block[i];
   }
   while (ia <= n_inc && ((H0 <= (t << 1)) || sum + (inc0.block[ia-1] << 1) + H0 <= 0))
   {
      dH = dH - inc0.block[ia - 1];
      if (dH > 0)t = dH;
      else t = -dH;
      ia++;
      sum += inc0.block[ia - 1 - 1];
   }



   while (ia < inc0.num - 1)
   {
      if (inc0.block[ia - 1] / 2 == inc0.block[ia] && abs(dH)>abs(inc0.block[ia - 1] + inc0.block[ia]))
      {
         inc[ia - 1] = inc0.block[ia - 1] * 2;
         dH = dH - inc[ia - 1];
         dh = find_basic_block(abs(dH));
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