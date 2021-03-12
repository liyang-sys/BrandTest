#include "modify_Hilbert3D.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

int **modify_Hilbert3D_ScanByBlock(int H, int W, block0 h_block, block0 w_block, int rank, int *f, int Hn0)
{
	int i, j;
	int Hn = H, Wn = w_block.block[iw] ;      
	int k0 = rank;
	int k = pow(2, rank);
	int g=0;
	int finish = f[0], lenw = w_block.num, lenh = h_block.num;
	int c_offset = 0;
	int idx;
	int **m0;
	int **m1 = create_arr0(3, H*Hn0*k);
	int sum=0, sum1=0;
	int dWn;
	int total=0;
	W1 = 0;
	Hn_next = Hn0, Wn_next = h_block.block[ih];
	for (i = ih; i<h_block.num; i++)
		sum = sum + h_block.block[i];
	total = H*w_block.block[iw] * k;
	while (Wn * 2 > Hn && (Wn_next < Hn_next || Wn_next == sum) && finish == 0){
		idx = w_block.order[iw];
		W1 = W1 + w_block.block[iw];
		if (idx>k0 || (idx == k0 && Hn <= k)){			
			m0 = modifyHilbert3D_BlockDeformation0(rank, idx, H);  
			for (i = total - H*w_block.block[iw] * k, j = 0; i < total; j++, i++){
				m1[0][i] = m0[0][j];
				m1[1][i] = m0[1][j] + c_offset;
				m1[2][i] = m0[2][j];
			}

			c_offset = c_offset + w_block.block[iw];  
			Hn_next = Hn_next - w_block.block[iw];
			if (iw+1 <w_block.num){
				iw = iw + 1;
				Wn = w_block.block[iw];
			}
			else{
				swap(&f[0], &f[1]);
				finish = f[0];
			}
		}
		else{
			for (i=iw; i<w_block.num; i++)   
				sum1 = sum1 + w_block.block[i];
			dWn = sum1;
			total = total + H * (dWn-w_block.block[iw] ) * k;
			m0 = zigzag3D(Hn, dWn, rank);
			for (i = iw + 1; i < w_block.num; i++)
				W1 = W1 + w_block.block[i];
			for (i = total-Hn* dWn*k,j=0; i < total ; j++,i++){
				m1[0][i] = m0[0][j];
				m1[1][i] = m0[1][j] + c_offset;
				m1[2][i] = m0[2][j];
			}
			swap(&f[0], &f[1]);
			finish = f[0];
		}
		total =total+ H*w_block.block[iw] * k;
		for (i = 0; i < 3; i++)
			free(m0[i]);
		free(m0);
	}
	total = total - H*w_block.block[iw] * k;
	total_block = total;
	return m1;
}

int *modify_Hilbert3D(int W, int H, int rank)
{
	int i, ia=0;
	int trsz = 0, tem, tem0;
	int k = pow(2, rank);
	block0 h_block, w_block;
	int Hn0;
	int r_offset = 0, c_offset = 0, mmax;
	int Hn , Wn;
	int **m, **m0, *m1;
	int vertical = 1;
	int f[2] = { 0, 1 };
	int finish = f[0];
	int tem1,  tem4;
	int *t, *tem2 = row0(20), *tem3 = row0(20), *tem5;
	int r_max = 0, c_max=0;
	int Hc, Wc;                           
	int total;
	int idx;
	m1 = (int*)calloc(((W*H) << rank)-1, sizeof(int));
	if (!m1){
		printf("创建m1数组失败！\n");
		exit(1);
	}
	ih = 0, iw = 0;
	if (H > W){
		tem = H; H = W; W = tem; trsz = 1;
	}
	Hn = H;
	Hn0 = W;
	Hc = H; Wc = W;
	h_block = find_basic_block0(H);
	w_block = find_basic_block0(W);
	Wn = w_block.block[iw];
 	m = modify_Hilbert3D_ScanByBlock(H, W, h_block, w_block,rank,f,Hn0);
	finish = f[0];
	m0 = create_arr0(3, H*W*k);
	Hn0 = Hn;
	mmax = m[1][0];
	for (i = 0; i < total_block; i++){
		m0[0][i] = m[0][i];
		m0[1][i] = m[1][i];
		m0[2][i] = m[2][i];
		if (mmax < m[1][i])
			mmax = m[1][i];
	}
	Hn0 = Hn;
	c_offset = c_offset + mmax;
	total = total_block;
	while (finish == 0){
		tem1 = H;
		for (i = 0; i < h_block.num; i++){
			tem2[i] = h_block.order[i];
			tem3[i] = h_block.block[i];
		}
		tem4 = ih;
		H = W;
		for (i = 0; i < w_block.num; i++){
			h_block.order[i] = w_block.order[i];
			h_block.block[i] = w_block.block[i];
		}
		ih = iw;
		W = tem1;
		for (i = 0; i < h_block.num; i++){
			w_block.order[i] = tem2[i];
			w_block.block[i] = tem3[i];
		}
		iw = tem4;
		tem0 = w_block.num;
		w_block.num = h_block.num;
		h_block.num = tem0;
		if (i >= w_block.num){
			w_block.order[i] = 0;
			w_block.block[i] = 0;
		}
		if (i >= h_block.num){
			h_block.order[i] = 0;
			h_block.block[i] = 0;
		}
		Hn = Hn_next; Wn = Wn_next;
		if (Hn >= 1){
			if (Wn >= k){
            if (m != NULL){                //注意释放前一次迭代的m的内存
               for (i = 0; i < 3; i++)
                  free(m[i]);
               free(m);
            }
				f[0] = 0; f[1] = 1;
				m = modify_Hilbert3D_ScanByBlock(Hn, Wn, h_block, w_block, rank, f,Hn0);
				finish = f[0];
			}
			else{
				Wn = 0;
				for (i = iw; i<w_block.num; i++)   
					Wn = Wn + w_block.block[i];
				W1 = Wn;

            if (m != NULL){
               for (i = 0; i < 3; i++)
                  free(m[i]);
               free(m);
            }
				m = zigzag3D(Hn, Wn, rank);
				swap(&f[0], &f[1]);
				finish = f[0];
			}
		}
		else{
			swap(&f[0], &f[1]);
			finish = f[0];
		}
		tem5 = row0(H*W*k);
		if (vertical > 0){
			idx = w_block.order[iw - 1];
			for (i = 0; i < W1*Hn*k; i++)
				tem5[i] = m[1][i];
			for (i = 0; i < W1*Hn*k; i++)
				m[1][i] = m[0][i]+c_offset;
			for (i = 0; i < W1*Hn*k; i++)
				m[0][i] = tem5[i]+r_offset;
			for (i = total; i < W1*Hn*k + total; i++)
				m0[0][i] = m[0][i - total];
			for (i = total; i < W1*Hn*k + total; i++)
				m0[1][i] = m[1][i - total];
			for (i = total; i < W1*Hn*k + total; i++)
				m0[2][i] = m[2][i - total];
			for (i = 0; i < W1*Hn*k; i++)			
			if (r_max < m[0][i]) r_max = m[0][i];
			r_offset = r_max;
			total = total + W1*Hn*k;
		}
		else{
			for (i = 0; i < W1*Hn*k; i++)
				m[1][i] = m[1][i] + c_offset;
			for (i = 0; i < W1*Hn*k; i++)
				m[0][i] = m[0][i] + r_offset;
			for (i = total; i < W1*Hn*k + total; i++)
				m0[0][i] = m[0][i - total];
			for (i = total; i < W1*Hn*k + total; i++)
				m0[1][i] = m[1][i - total];
			for (i = total; i < W1*Hn*k + total; i++)
				m0[2][i] = m[2][i - total];
			for (i = 0; i < W1*Hn*k; i++)                           
			if (c_max < m[1][i])
				c_max = m[1][i];
				c_offset = c_max;
				total = total + W1*Hn*k;
				W1 = Hn;
		}
		Hn0 = Hn;
		vertical = -vertical;
		free(tem5);
	}
   if (m != NULL){
      for (i = 0; i < 3; i++)
         free(m[i]);
      free(m);
   }
	if (trsz == 1){
		t = row0(Hc*Wc*k);
		for (i = 0; i < Hc*Wc*k; i++){
			t[i] = m0[1][i];
			m0[1][i] = m0[0][i];
			m0[0][i] = t[i];
		}
		free(t);
	}
   
	for (i = 0; i < (((H*W)<<rank) - 1); i++){   /*将二维的m变为一维的m1*/
		if (m0[2][i + 1] - m0[2][i] == 1)/*时间维*/
			m1[i] = 1;
		else if (m0[2][i + 1] - m0[2][i] == -1)
			m1[i] = -1;
		else if (m0[1][i + 1] - m0[1][i] == 1)/*列*/
			m1[i] = 3;
		else if (m0[1][i + 1] - m0[1][i] == -1)
			m1[i] = -3;
		else if (m0[0][i + 1] - m0[0][i] == 1)/*行*/
			m1[i] = 2;
		else
			m1[i] = -2;
	}

	free(m0[0]); free(m0[1]); free(m0[2]); free(m0);
	free(tem2); free(tem3);             //释放内存
	return m1;
}



