#include "modify_Hilbert2D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//该函数用于实现向下加的功能，传入已有数组以及它的长和宽，以及要添加的长度，输出添加完毕的数组
int **add(int **arr, int width, int height, int num)
{
	int i, j, k;
	int counter = height*width;
	int flag = 0;
	int offset;
	int portion;
	int idy, idx,pt,pt_new;
	int unit = num * 2;
	int new_unit = num * 2 + num;
	int orginal = height*width;
	int sum = width*(height+num);
	int pointer, pointer0;
	int subnum = width / (num*2);

	int *st = (int *)malloc((subnum+1)*sizeof(int *));
	int *orginal_row = (int *)malloc(orginal*sizeof(int *));								//用于存储传入数组的横坐标
	int *orginal_culumn = (int *)malloc(orginal*sizeof(int *));								//用于存储传入数组的纵坐标
	int *sum_row = (int *)malloc(sum*sizeof(int *));										//用于存储完成加法后的数组的横坐标
	int *sum_culumn = (int *)malloc(sum*sizeof(int *));										//用于存储完成加法后的数组的横坐标
	int *transit_row = (int *)malloc(sum*sizeof(int *));									//过渡数组,与sum长度相等
	int *transit_culumn = (int *)malloc(sum*sizeof(int *));									//过渡数组,与sum长度相等
		
	int *row = (int *)malloc((unit*unit)*sizeof(int *));									//用于存储所添加模块的横坐标
	int *culumn = (int *)malloc((unit*unit)*sizeof(int *));									//用于存储所添加模块的纵坐标
	int *new_row = (int *)malloc((new_unit*unit)*sizeof(int *));							//用于存储添加后模块的横坐标
	int *new_culumn = (int *)malloc((new_unit*unit)*sizeof(int *));							//用于存储添加后模块的纵坐标

	int **orginal_arr = (int **)malloc(width*sizeof(int *));								//用于表示输入数组
	int **new_arr = (int **)malloc(width*sizeof(int *));									//用于表示输出数组
	for (int i = 0; i<width; i++)
		new_arr[i] = (int *)malloc((height+num)*sizeof(int));								
	for (int i = 0; i<width; i++)
		orginal_arr[i] = (int *)malloc(height*sizeof(int));									
///////////////////////////////////////////////////////////////////////////////////////////////初始定义结束
	//*****如果所增加高度与原高度相等，则采用另一种加法
	if (width == height && num == height){													
		idx = 0;
		idy = 0;
		offset = num*num / 4;
		portion = num*num / 4;
		for (j = 0; j < height/2; j++){														//第一阶段，与原第一部分一致
			for (i = 0; i < width/2; i++)
				new_arr[i][j] = arr[i][j];
		}
		idy = height / 2;
		for (j = height/2; j < height; j++){												//第二阶段，与原第一部分一致加上offset
			for (i = 0; i < width / 2; i++)
				new_arr[i][j] = arr[i][j - idy] + offset;
		}
		idy = height;
		offset = offset + portion;
		for (j = height; j < height + num; j++){											//第三阶段，与原数组整个一直
			for (i = 0; i < width; i++)
				new_arr[i][j] = arr[i][j - idy] + offset;
		}
		offset = offset + portion;
		idy = height / 2;
		for (j = height / 2; j < height; j++){												//第四阶段，与原第四部分一致加上offset
			for (i = width / 2; i < width; i++)
				new_arr[i][j] = arr[i][j - idy] + offset;
		}
		offset = offset + portion;
		for (j = 0; j < height / 2; j++){													//第五阶段，与原第四部分一致加上offset
			for (i = width / 2; i < width; i++)
				new_arr[i][j] = arr[i][j] + offset;
		}
	}
	//********如果需要增加的部分比原高度小，则使用常规加法
	else{
		//用于得到传进来有用部分数组
		for (i = 0; i < width; i++){
			for (j = 0; j < height; j++)
				orginal_arr[i][j] = arr[i][j];
		}
		idy = 0;
		idx = 0;
		for (i = 0; i < unit*unit; i++)															//这是找到要用来添加的单元
			d2xy(unit, i, &row[i], &culumn[i]);
		//////////////开始拼起来一个新的
		offset = unit*unit / 4;
		portion = unit*unit / 4;
		pt = 0;
		pt_new = 0;
		for (i = pt_new; i < offset; i++){														//第一阶段，与原来一致
			new_row[i] = row[i];
			new_culumn[i] = culumn[i];
		}
		idx += 0;
		idy += num;
		pt_new = offset;
		offset += portion;
		for (i = pt_new; i < offset; i++){														//第二阶段，与上一阶段一样，x相同，y+num
			new_row[i] = new_row[i - num*num];
			new_culumn[i] = new_culumn[i - num*num] + num;
		}
		pt_new = offset;
		pt += portion;
		offset += 2 * portion;
		idx += 0;
		idy += num;
		for (i = pt_new, j = pt; i < offset; i++, j++){											//第三阶段，等于原单元的第二部分
			new_row[i] = row[j];
			new_culumn[i] = culumn[j] + num;
		}
		pt_new = offset;
		pt += 2 * portion;
		idx += num;
		idy -= num;
		offset += portion;
		for (i = pt_new, j = pt; i < offset; i++, j++) {										//第四阶段，与原单元格第三部分一致，y+num;
			d2xy(unit, j, &new_row[i], &new_culumn[i]);
			new_culumn[i] += num;
		}
		pt_new = offset;
		idx += 0;
		idy -= num;
		offset += portion;
		for (i = pt_new, j = pt; i < offset; i++, j++) {										//第五阶段，与原单元格第三部分一致
			new_row[i] = row[j];
			new_culumn[i] = culumn[j];
		}
		/////////需要计算出每个要加的结点的位置在哪
		idx = 0;
		idy = height - num * 2;
		for (i = 0; i < subnum; i++, idx += unit)
			st[i] = orginal_arr[idx][idy];
		//此处是把传进来的那个二维数组转化成行跟列的两个对应x和y坐标的一维数组
		for (i = 0; i < width; i++){
			for (j = 0; j < height; j++){
				orginal_culumn[orginal_arr[i][j]] = j;
				orginal_row[orginal_arr[i][j]] = i;
			}
		}
		//使过渡函数等于原始函数
		for (i = 0; i < width*height; i++){
			transit_row[i] = orginal_row[i];
			transit_culumn[i] = orginal_culumn[i];
		}
		idx = 0;
		idy = height - num * 2;
		/////////接下来要把其他的全部改变一次，条件是单元格自加等于len
		for (i = 0; i < subnum; i++){															//循环添加	
			pointer = st[i];
			st[i + 1] = st[i + 1] + unit*num*(i + 1);
			counter += num*unit;
			for (j = 0; j < pointer; j++){														//第一部分与原数组一致
				sum_row[j] = transit_row[j];
				sum_culumn[j] = transit_culumn[j];
			}
			pointer += new_unit*unit;
			for (j = st[i], k = 0; j < pointer; j++, k++){										//第二部分把已修改部分添加进去
				sum_row[j] = new_row[k] + idx;
				sum_culumn[j] = new_culumn[k] + idy;
			}
			pointer0 = st[i] + unit*unit;
			for (j = pointer, k = pointer0; j < counter; j++, k++){								//第三部分与原数组一致
				sum_row[j] = transit_row[k];
				sum_culumn[j] = transit_culumn[k];
			}
			idx += unit;
			for (j = 0; j < sum; j++){															//使过渡数组等于当前数组用于下一次计算
				transit_row[j] = sum_row[j];
				transit_culumn[j] = sum_culumn[j];
			}
		}
		//把两个一维数组对应的x和y的坐标恢复成二维数组
		for (k = 0; k < width*(height + num); k++)
			new_arr[sum_row[k]][sum_culumn[k]] = k;
	}
	return new_arr;
}

//该函数用于实现往上减的功能，传入二维数组以及它的长和宽，以及要减去的长度
int **sub(int **arr, int width, int height, int num)
{
	int i, j, k;
	int unit = num * 2;
	int unit_size = unit*unit;
	int new_unit = num * 2 - num;
	int new_unit_size = new_unit*unit;
	int subnum = width / unit;
	int pt, pt_new, pointer, pointer0;
	int counter = width*height;
	int idx, idy;
	int orginal = width*height;
	int sum = width*(height - num);

	int *st = (int *)malloc((subnum + 1)*sizeof(int *));
	int *orginal_row = (int *)malloc(orginal*sizeof(int *));									//原来的总长哦
	int *orginal_culumn = (int *)malloc(orginal*sizeof(int *));
	int *transit_row = (int *)malloc(orginal*sizeof(int *));									//过渡数组,与orginal一样长
	int *transit_culumn = (int *)malloc(orginal*sizeof(int *));
	int *transit_row1 = (int *)malloc(orginal*sizeof(int *));									//过渡数组,与orginal一样长
	int *transit_culumn1 = (int *)malloc(orginal*sizeof(int *));
	int *sum_row = (int *)malloc(sum*sizeof(int *));											//创建总的新的数组的总长的数组
	int *sum_culumn = (int *)malloc(sum*sizeof(int *));

	int *row = (int *)malloc((unit_size)*sizeof(int *));										//创建要sub的那个unit的x坐标的数组？
	int *culumn = (int *)malloc((unit_size)*sizeof(int *));										//创建要sub的那个unit的y坐标的数组？
	int *new_row = (int *)malloc((new_unit_size)*sizeof(int *));								//创建原来的那个unit并加上部分的x坐标的数组？
	int *new_culumn = (int *)malloc((new_unit_size)*sizeof(int *));								//创建原来那个unit加上部分的y坐标的数组？

	int **orginal_arr = (int **)malloc(width*sizeof(int *));									//创建行
	for (int i = 0; i<width; i++)
		orginal_arr[i] = (int *)malloc(height*sizeof(int));										//创建列

	//用于得到传进来有用部分数组
	for (i = 0; i < width; i++){
		for (j = 0; j < height; j++)
			orginal_arr[i][j] = arr[i][j];
	}
	for (i = 0; i < unit_size; i++)																//这是找到要用来添加的单元
		d2xy(unit, i, &row[i], &culumn[i]);
	idx = 0;
	idy = -1*(num);
	//开始拼起来一个新的，把1和3部分删除
	pt = unit_size / 4; 
	pt_new = pt*3;
	for (i = pt, j = 0; i < pt_new; i++, j++){
		new_row[j] = row[i] + idx;
		new_culumn[j] = culumn[i] + idy;
	}
	//先把要进行计算的位置找出来
	idx = 0;
	idy = height - num * 2;
	for (i = 0; i < subnum; i++, idx += unit)
		st[i] = orginal_arr[idx][idy];
	//此处是把传进来的那个二维数组转化成行跟列的两个对应x和y坐标的一维数组
	for (i = 0; i < width; i++){
		for (j = 0; j < height; j++){
			orginal_culumn[orginal_arr[i][j]] = j;
			orginal_row[orginal_arr[i][j]] = i;
		}
	}
	for (i = 0; i < width*height; i++){															//初始化过渡数组，令其与原数组相同
		transit_row[i] = orginal_row[i];
		transit_culumn[i] = orginal_culumn[i];
		transit_row1[i] = orginal_row[i];
		transit_culumn1[i] = orginal_culumn[i];
	}
	//开始把新单元装回到原来的大数组里
	idx = 0;
	idy = height - unit;
	/////////接下来要把其他的全部改变一次，条件是单元格自加等于len
	for (i = 0; i < subnum; i++){																//循环减
		pointer = st[i];
		st[i + 1] = st[i + 1] - new_unit_size*(i + 1);
		counter -= new_unit_size;
		for (j = 0; j < pointer; j++){															//与原数组第一部分一致
			transit_row1[j] = transit_row[j];
			transit_culumn1[j] = transit_culumn[j];
		}
		pointer += new_unit_size;
		for (j = st[i], k = 0; j < pointer; j++, k++){											//第二部分用已修改单元覆盖
			transit_row1[j] = new_row[k] + idx;
			transit_culumn1[j] = new_culumn[k] + idy;
		}
		pointer0 = st[i] + unit_size;
		for (j = pointer, k = pointer0; j < counter; j++, k++){									//第三部分与原第三部分一致
			transit_row1[j] = transit_row[k];
			transit_culumn1[j] = transit_culumn[k];
		}
		idx += unit;
		for (j = 0; j < orginal; j++){																
			transit_row[j] = transit_row1[j];
			transit_culumn[j] = transit_culumn1[j];
		}
	}
	for (i = 0; i < sum; i++){																	//获得最后的数组
		sum_row[i] = transit_row1[i];
		sum_culumn[i] = transit_culumn1[i];
	}
	//把两个一维数组对应的x和y的坐标恢复成二维数组
	for (k = 0; k < sum; k++)
		arr[sum_row[k]][sum_culumn[k]] = k;
	return arr;
}

