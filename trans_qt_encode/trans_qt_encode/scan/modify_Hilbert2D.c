/*
*修改日期：2021-1-13
*bug：1,calculation 函数计算分块时，对减法操作划分出的子块个数计算错误
*     2,judge_order 重新计算分块时，对新分块个数以及阶数计算错误
*
*修改日期：2021-1-14
*改善：路径用行和列的坐标表示，简化了转换过程
*
*修改日期：2021-1-15---2021-1-18
*逻辑上的修改：重新划分子块，重新写lfind_2s_inc函数，就是find_2s_inc_2
*              使整个扫描结果和MATLAB中一致
*              对draw_hibert和calculation都有相应改动
*
*----张霞
*/


#include "modify_Hilbert2D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*用于动态创建任意长度二维数组并返回*/
int **create_arr(int width, int height)
{
	int **new_arr = (int **)malloc(width*sizeof(int *));										//创建行
	for (int i = 0; i<width; i++)
		new_arr[i] = (int *)malloc(height*sizeof(int));											//创建列
	return new_arr;
}

//传递基本square的长度width,和最后要加成的高度height，还有h要加的序列，返回一个加完的二维数组
int **bc_operation(int width, int height, block h)
{
	int i, j;
	int len = width;
	int times = h.num;
	int **block_one = create_arr(len, height);
	int **square = create_arr(len, len);
	int idh = width;
	int w = len; 
	for (i = 0; i < len; i++){
		for (j = 0; j < len; j++)
			square[i][j] = xy_point2d(len, i, j);
	}

	for (i = 0; i < times; i++){
		if (h.block[i] > 0){
			square = add(square, len, idh, h.block[i]);
			idh += h.block[i];
		}
		else if (h.block[i] < 0){
			square = sub(square, len, idh, abs(h.block[i]));
			idh += h.block[i];
		}
	}
	for (i = 0; i < height; i++){
		for (j = 0; j < len; j++)
			block_one[j][i] = square[j][i];
	}
	return block_one;
}


total calculation(int w, int h)
{
	total basic;
   int MaxDel = 0, ia = 0;
   block th_block;

	basic.basic_block = find_basic_block(w);													//计算出关于宽度的序列
	basic.h_order = find_basic_block(h);														//计算出关于高度的序列
   basic = judge_order_2(basic, w, h);
	//basic = judge_order(basic, w, h);															//结合高度和宽度的序列判断并修改已有序列
	//basic.h_order = find_2s_inc(basic.h_order, h - basic.basic_block.block[0]);					//计算高度的序列是否能够缩短
      
   th_block = find_2s_inc_2(basic.basic_block.block[0], h - basic.basic_block.block[0], basic.h_order);
   
   while (ia<th_block.num&&th_block.block[ia]>0){
      MaxDel += th_block.block[ia];
      ia++;
   }
   if (basic.basic_block.block[0] == MaxDel){
      return basic;
   }
   else{
      int i = 0;
      while (th_block.block[i] != 0){
         basic.h_order.block[i] = th_block.block[i];
         i++;
      }
      while (basic.h_order.block[i] != 0){
         basic.h_order.block[i] = 0;
         i++;
      }
      basic.h_order.num = th_block.num;
   }

	return basic;
}

int **draw_hibert(int W, int H)
{
	int flip = 0;;
	int i, j, k;
	total chunk;
   block h_new, tem_h1, tem_w1;
	int w , h;																					//用于表示每次要计算的w和h
	int W1, H1, W2, H2;
	int w1, h1, w2, h2;																			//高度和宽度调转时使用
	int flag = 0;																				//用于标志是否遍历，只会变一次
   int ptx = 0, pty = 0, ptw = 0, ptz = 0;																//用于作为标记要添加该数组的位置指针
	int offset = 0;																				//表示此时已添加的数
   int rotation = 0, flg1 = 1, flg2 = 1, f1 = 0, f2 = 0;			//用于判断数组是否需要旋转，0表示旋转0度，1表示顺时针旋转270度,flg一起判断是否改变方向,f辅助计算坐标加的次数
	int **img;																					//这个是最后的总的图像数组
	int **part, **part_rot, **part_culumn;														//每次用来计算的新数组
	int *part_row;			
	int **img_flip;

	if (W < H){																					//如果初始宽度小于高度，把宽和高调换处理
		flip = 1;
		w = H;
		h = W;
		W2 = H;
		H2 = W;
	}
	else{
		flip = 0;
		w = W;
		h = H;
		W1 = W;
		H1 = H;
	}
	img = create_arr(w, h);
	w1 = w;
	h1 = h;
	while (flag == 0){
		if ((rotation == 0||flg1==0)&&flg2!=0){																		//先判断是否需要翻转
         f2 = 0;
         if (f1 > 0)pty -= h1;
         if (w1 < h1&&flg1==1){
				rotation = 1;
				h2 = w1;
				w2 = h1;
			}
			else if (w1 >= h1 && h1 == 1){														//判断是否会遇到高度只有一行的情况
				part_row = row(w1);
				for (i = 0; i < w1; i++)
					img[ptx + i][pty] = part_row[i] + offset;
				ptx += w1;
				pty += 1;
            ptz = pty;
			}
			else if ( h1 > 1){										//w1 >= h1 &&	//高度大于1
				chunk = calculation(w1, h1);	
				for (i = 0; i < chunk.basic_block.num; i++){
					part = bc_operation(chunk.basic_block.block[i], h1, chunk.h_order);
					for (j = 0; j < chunk.basic_block.block[i]; j++){
						for (k = 0; k < h1; k++)
							img[ptx + j][pty + k] = part[j][k] + offset;
					}
               ptz = pty + h1;

					for (int ii = 0; ii < chunk.basic_block.block[i]; ii++)					
						free(part[ii]);														
					free(part);																	
					offset += chunk.basic_block.block[i] * h1;
					ptx += chunk.basic_block.block[i];               
					if (chunk.basic_block.block[i + 1] != chunk.basic_block.block[i]){
						i += 1;
						break;
					}
				}
            flg1 = 1;
				rotation = 1;																	//方向顺时针改变270度
				w -= i*chunk.basic_block.block[0];												//现在的宽度变化
				w2 = h1;																		//下一个过程的宽度
				h2 = w;																			//下一个过程的高度
            tem_h1 = find_basic_block(h1);
            if ((chunk.basic_block.block[0] << 1) > h1&&tem_h1.block[0] <= w){
               flg1 = 0;
               w1 = w;
               f1++;
               pty += h1;
               ptz = pty;
            }

			}
			if (flip == 0){
				if (ptx == W && ptz == H)														//用于判断是否遍历完毕
					flag = 1;
				else
					flag = 0;
			}
			else{
				if (ptx == H && ptz == W)														//用于判断是否遍历完毕
					flag = 1;
				else
					flag = 0;
			}
		}
		else {    //if (rotation == 1)
         if (f2 > 0)ptx -= w1;;
         f1 = 0;
			if (h2 == 1){
				part_culumn = create_arr(h2, w2);
				for (i = 0; i < w2; i++)
					part_culumn[0][i] = i;
				for (i = 0; i < w2; i++)
					img[ptx][pty + i] = part_culumn[0][i] + offset;
				pty += w2;
				ptx += 1;
				ptw = ptx;
			}
			else if (h2 > 1){
				h_new = find_basic_block(w2);
				w2 = h_new.block[0];
				chunk = calculation(w2, h2);
				for (i = 0; i < chunk.basic_block.num; i++){
					part = bc_operation(chunk.basic_block.block[i], h2, chunk.h_order);
					part_rot = rotate_arr(chunk.basic_block.block[i], h2, part);
					for (j = 0; j < chunk.basic_block.block[i]; j++){
						for (k = 0; k < h2; k++)
							img[ptx + k][pty + j] = part_rot[k][j] + offset;
					}
					ptw = ptx + h2;
					//释放两个二维数组
					for (int ii = 0; ii < chunk.basic_block.block[i]; ii++)						//释放内存
						free(part[ii]);															//释放列
					free(part);																	//释放行
					for (int ii = 0; ii < h2; ii++)												//释放内存
						free(part_rot[ii]);														//释放列
					free(part_rot);																//释放行

					offset += chunk.basic_block.block[i] * h2;
					pty += chunk.basic_block.block[i];
					if (chunk.basic_block.block[i + 1] != chunk.basic_block.block[i])
						break;
				}
				rotation = 0;																	//方向顺时针改变90度
            flg2 = 1;
            h -= w2;
				w1 = h2;																		//下一个宽度
				h1 = h;																			//下一个高度
            tem_h1 = find_basic_block(h1);
            tem_w1 = find_basic_block(w1);
            if ((tem_h1.block[0] << 1) > h2&&tem_w1.block[0] <= h1){
               flg2 = 0;
               w2 = h;
               f2++;
               ptx += w1;
               ptw = ptx;
            }

			}
			else if (h2 == 0){
				flag = 1;
				break;
			}
			if (flip == 0){
				if (ptw == W && pty == H)														//用于判断是否遍历完毕
					flag = 1;
				else
					flag = 0;
			}
			else if (flip == 1){
				if (ptw == H && pty == W)														//用于判断是否遍历完毕
					flag = 1;
				else
					flag = 0;
			}
		}      

	}
	if (flip == 1){
		img_flip = rotate_arr(H, W, img);
		return img_flip;
	}
	else if (flip == 0)
		return img;
}

//生成一维最后一行数组
int *row(int n)
{
	int *arr = (int *)malloc(n*sizeof(int *));
	for (int i = 0; i < n; i++)
		arr[i] = i;
	return arr;
}


//////用于旋转二维数组，w为原数组宽，h为原数组高
int **rotate_arr(int w, int h, int **old)
{
	int i, j;
	int **new, **mirror;
	new = create_arr(h, w);
	mirror = HorizontalMirror(old, w, h);				//镜像对称后再旋转
	////////270,old[w][h],new[h][w]
	for (i = 0; i<h; i++)
	{
		for (j = 0; j<w; j++)
			new[i][j] = mirror[w - 1 - j][i];
	}
	return new;
}

//二维数组水平镜像
int **HorizontalMirror(int **arr, int n, int m)
{
	int r, c;
	for (r = 0; r < m; r++){
		for (c = 0; c < n / 2; c++)						//每一行转换（n/2）次
		{												//将左右两个对称的元素进行交换	
			int temp = arr[c][r];
			arr[c][r] = arr[n-c-1][r];
			arr[n-c-1][r] = temp;
		}
	}
	return arr;
}

int *modify_Hilbert2D(int W, int H)
{
	int i, j, k, lg;
	int **m0, **img, *m;
	lg = H*W;

	img = draw_hibert(W, H);

   m0 = create_arr(3, lg);     //把路径用行和列的坐标表示;m0[0]:行,m0[1]:列,m0[3]:深度,m0[3]在二维中全为0；
   for (i = 0; i < lg; i++)
      m0[2][i] = 0;

   int x = 0, y = 0;
   m0[0][0] = 0, m0[1][0] = 0;
   for (int i0 = 1; i0 < lg; i0++){
      if (x>0 && img[x - 1][y] == img[x][y] + 1){
         x = x - 1;
         m0[0][i0] = y;
         m0[1][i0] = x;
      }
      else if (x<W - 1 && img[x + 1][y] == img[x][y] + 1){
         x = x + 1;
         m0[0][i0] = y;
         m0[1][i0] = x;
      }
      else if (y>0 && img[x][y - 1] == img[x][y] + 1){
         y = y - 1;
         m0[0][i0] = y;
         m0[1][i0] = x;
      }
      else if (y<H - 1 && img[x][y + 1] == img[x][y] + 1){
         y = y + 1;
         m0[0][i0] = y;
         m0[1][i0] = x;
      }

   }
   for (i = 0; i < W; i++)
      free(img[i]);
   free(img);

   /*************************/
	//for (int i = 5524; i < 5534; i++)
	//	printf("%d ", m0[0][i]);
	//printf("\n");
	//for (int i = 5524; i < 5534; i++)
	//	printf("%d ", m0[1][i]);

	//FILE *fp2;
	//fopen_s(&fp2, "m0.dat", "wb");
	//for (i = 0; i<2; i++)
	//{
	//	fwrite(&m0[i][0], sizeof(int), lg, fp2);
	//}
	//fclose(fp2);
	/*************/

   m = (int*)calloc(lg, sizeof(int));
   if (!m){
      printf("创建数组m失败！\n");
      exit(1);
   }
	for (i = 0; i < lg - 1; i++){   /*将二维的m变为一维的m1*/
		if (m0[0][i + 1] - m0[0][i] == 1)/*行*/
			m[i] = 2;
		else if (m0[0][i + 1] - m0[0][i] == -1)
			m[i] = -2;
		else if (m0[1][i + 1] - m0[1][i] == 1)/*列*/
			m[i] = 3;
		else if (m0[1][i + 1] - m0[1][i] == -1)
			m[i] = -3;
	}
	free(m0[0]); free(m0[1]); free(m0[2]); free(m0);

	return m;
}

