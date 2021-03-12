#include "modify_Hilbert3D.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

//传入二维数组和已经构造好的路径，以及要添加的长度
int **modify_Hilbert3D_sub(int dim, int inc, int **m0, int rank, int *st0)
{                                                                                 
	int i, j, ii, jj, k = dim*dim*pow(2, rank), k2 = dim*(dim + inc)*pow(2, rank);
	int  pt, pt0; 
	int subdim = abs(inc) * 2, nsub = dim / subdim,  L1;   
	int *st, *stnew;   
	int hN, hN1, qN, qN1, hqN, hqN1;
	int Nsub1, Nq, Nq1, Nh;
	int *rsub, *csub, *dsub;
	int len;           
	int ed0 = 1, ednew0 = 1, s,ed;
	int *rsub0, *csub0, *dsub0, *sr3, *sc3, *sd3;   
	int *sr1, *sc1, *sd1, *sr2;
	int *sr5, *sc5, *sd5, *sr4;
	int snew, ednew;
	int *sc2, *sd2, *sc4, *sd4;
	int **m;
	int L0= L;
	int max = m0[0][0], min = m0[0][0];
	for (int p = 0; p < 3; p++)
	for (j = 0; j < N0; j++)
	if (max < m0[p][j])max = m0[p][j];
	for (int p = 0; p < 3; p++)
	for (j = 0; j < N0; j++)
	if (min > m0[p][j])min = m0[p][j];
	if (dim < subdim)
		nsub++;
	st = row0(nsub); stnew = row0(nsub);
	for (i = 0; i < nsub; i++)
		st[i] = 1;
	L1 = log2(nsub);
	L = L1;
	for (i = L0; i <= L1; i++)                //起始点大于1的情况（找起始点）
	{
		nsubtem = nsubtem * 2;
		Nsub = Nsub / 4;
		for (pt = 1; pt < nsubtem; pt += 2)
		{
			pt0 = pt / 2;
			st[pt - 1] = st0[pt0] + Nsub;
			st[pt] = st[pt - 1] + Nsub;
		}
		for (j = 0; j < nsub; j++)
			st0[j] = st[j];
	}
	hN = Nsub / 2; hN1 = hN + 1; qN = Nsub / 4; qN1 = qN + 1; hqN = hN + qN; hqN1 = hqN + 1;
	Nsub1 = Nsub + 1; Nq = Nsub + qN; Nq1 = Nq + 1; Nh = Nsub + hN;
   if (inc>0){
      if (dim - inc > 0){
         len = N0 + nsub*hN;
         m = create_arr0(3, len);
         for (i = 0; i < nsub; i++){
            j = i*hN;
            stnew[i] = st[i] + j;
         }
      }
      else{
         len = N0 + N0;
         m = create_arr0(3, len);
         for (i = 0; i < nsub; i++){
            j = i*(hN + qN);
            stnew[i] = st[j] + j;
         }
      }
   }
	else{
		len = N0 - nsub*hN;
		m = create_arr0(3, len);
		for (i = 0; i < nsub; i++){
			j = i*hN;
			stnew[i] = st[i] - j;
		}
	}
	for (i = 0; i < nsub; i++)      //add
	{
		s = st[i]; ed = s + Nsub - 1;
		rsub0 = row0(ed - s + 1);
		csub0 = row0(ed - s + 1);
		dsub0 = row0(ed - s + 1);
		ii = 0;
		for (j = s-1; j < ed ; j++){
			rsub0[ii] = m0[0][j];
			csub0[ii] = m0[1][j];
			dsub0[ii] = m0[2][j];
			ii += 1;
		}
		sr3 = row0(hqN - qN1 + 1);
		sc3 = row0(hqN - qN1 + 1);
		sd3 = row0(hqN - qN1 + 1);
		ii = 0;
		for (j = qN1 - 1; j < hqN; j++){
			sr3[ii] = rsub0[j]+inc;
			sc3[ii] = csub0[j];
			sd3[ii] = dsub0[j];
			ii += 1;
		}
		if (inc>0){
			sr1 = row0(qN);
			sc1 = row0(qN);
			sd1 = row0(qN);
			sr5 = row0(Nsub - hqN1 + 1);
			sc5 = row0(Nsub - hqN1 + 1);
			sd5 = row0(Nsub - hqN1 + 1);
			if (dim - inc > 0)                  //There will be 6 sub-squares, and they are devided into 5 connected parts below.
			{
				rsub = row0(Nh);
				csub = row0(Nh);
				dsub = row0(Nh);
				sr2 = row0(qN);
				sr4 = row0(qN);
				for (j = 0; j < qN; j++){
					sr1[j] = rsub0[j];
					sc1[j] = csub0[j];
					sd1[j] = dsub0[j];
					sr2[j] = sr1[j] + inc;
				}
				ii = 0;
				for (j = hqN1 - 1; j < Nsub; j++){
					sr5[ii] = rsub0[j];
					sc5[ii] = csub0[j];
					sd5[ii] = dsub0[j];
					sr4[ii] = sr5[ii] + inc;
					ii += 1;
				}
				
				for (j = 0; j < Nh; j++){
					if (j < qN){
						rsub[j] = sr1[j];
						csub[j] = sc1[j];
						dsub[j] = sd1[j];
					}
					else if (j < hN && j >= qN1 - 1){
						rsub[j] = sr2[j - qN1 + 1];
						csub[j] = sc1[j - qN1 + 1];
						dsub[j] = sd1[j - qN1 + 1];

					}
					else if (j < Nsub && j >= hN1 - 1){
						rsub[j] = sr3[j - hN1 + 1];
						csub[j] = sc3[j - hN1 + 1];
						dsub[j] = sd3[j - hN1 + 1];
					}
					else if (j < Nq && j >= Nsub1 - 1){
						rsub[j] = sr4[j - Nsub1 + 1];
						csub[j] = sc5[j - Nsub1 + 1];
						dsub[j] = sd5[j - Nsub1 + 1];
					}
					else if (j < Nh && j >= Nq1 - 1){
						rsub[j] = sr5[j - Nq1 + 1];
						csub[j] = sc5[j - Nq1 + 1];
						dsub[j] = sd5[j - Nq1 + 1];
					}
				}
				snew = stnew[i]; ednew = stnew[i] + Nsub + hN - 1;
				ii = ed0-1;
				jj = 0;
				for (j = 0; j < ednew; j++){
					if (j >= ednew0 - 1 && j < snew - 1){
						m[0][j] = m0[0][ii];
						m[1][j] = m0[1][ii];
						m[2][j] = m0[2][ii];
						ii += 1;
					}
					else if (j >= snew - 1 && j < ednew){
						m[0][j] = rsub[jj];
						m[1][j] = csub[jj];
						m[2][j] = dsub[jj];
						jj += 1;
					}
				}
				ed0 = ed;  ednew0 = ednew;
				free(sr2);  free(sr4);free(rsub); free(csub); free(dsub);
			}
			else{
				rsub = row0(Nsub * 2);
				csub = row0(Nsub * 2);
				dsub = row0(Nsub * 2);
				sr2 = row0(2 * qN);
				sc2 = row0(2 * qN);
				sd2 = row0(2 * qN);
				sr4 = row0((Nsub - hqN1 + 1) * 2);
				sc4 = row0((Nsub - hqN1 + 1) * 2);
				sd4 = row0((Nsub - hqN1 + 1) * 2);
				for (j = 0; j < qN; j++){
					sr1[j] = rsub0[j];
					sc1[j] = csub0[j];
					sd1[j] = dsub0[j];
					sr2[j] = sr1[j] + inc / 2;
					sc2[j] = sc1[j];
					sd2[j] = sd1[j];
				}
				for (j = qN; j < 2 * qN; j++){
					sr2[j] = sr1[j - qN] + inc;
					sc2[j] = sc1[j - qN];
					sd2[j] = sd1[j - qN];
				}
				ii = 0;
				for (j = hqN1 - 1; j < Nsub; j++){
					
					sr5[ii] = rsub0[j];
					sc5[ii] = csub0[j];
					sd5[ii] = dsub0[j];
					sr4[ii] = sr5[ii] + inc;
					sc4[ii] = sc5[ii];
					sd4[ii] = sd5[ii];
					ii += 1;
				}
				for (j = ii; j < (Nsub - hqN1 + 1) * 2; j++){
					sr4[j] = sr4[j - ii ] - inc / 2;
					sc4[j] = sc5[j - ii ];
					sd4[j] = sd5[j - ii ];
				}
				for (j = 0; j < Nsub*2; j++){
					if (j < qN){
						rsub[j] = sr1[j];
						csub[j] = sc1[j];
						dsub[j] = sd1[j];
					}
					else if (j < qN + hN && j >= qN1 - 1){
						rsub[j] = sr2[j - qN1 + 1];
						csub[j] = sc2[j - qN1 + 1];
						dsub[j] = sd2[j - qN1 + 1];
					}

					else if (j < qN + Nsub && j >= hN1 + qN - 1){
						rsub[j] = sr3[j - hN1 - qN + 1];
						csub[j] = sc3[j - hN1 - qN + 1];
						dsub[j] = sd3[j - hN1 - qN + 1];
					}
					else if (j < Nq + hN && j >= Nsub1 + qN - 1){
						rsub[j] = sr4[j - Nsub1 - qN + 1];
						csub[j] = sc4[j - Nsub1 - qN + 1];
						dsub[j] = sd4[j - Nsub1 - qN + 1];
					}
					else if (j < Nsub * 2 && j >= Nq1+hN - 1){
						rsub[j] = sr5[j - Nq1 - hN + 1];
						csub[j] = sc5[j - Nq1 - hN + 1];
						dsub[j] = sd5[j - Nq1 - hN + 1];
					}
				}
				snew = stnew[i]; ednew = stnew[i] + Nsub + Nsub - 1;
				ii = ed0;
				jj = 0;
				for (j = 0; j < ednew; j++){
					if (j >= ednew0 - 1 && j < snew - 1){
						m[0][i] = m0[0][ii];
						m[1][i] = m0[1][ii];
						m[2][i] = m0[2][ii];
						ii += 1;
					}
					else if (j >= snew - 1 && j < ednew){
						m[0][j] = rsub[jj];
						m[1][j] = csub[jj];
						m[2][j] = dsub[jj];
						jj += 1;
					}
				}
				ed0 = ed;  ednew0 = ednew;
				free(sr2);free(sc2); free(sd2);  
				free(rsub); free(csub); free(dsub);
				free(sr4);free(sc4);free(sd4);
			}
			free(sr1); free(sc1); free(sd1);
			free(sr5); free(sc5); free(sd5);	
		}
		else{
			snew = stnew[i]; ednew = stnew[i] + hN - 1;
			ii = ed0 - 1;
			for (j = ednew0-1; j < snew - 1; j++){
				m[0][j] = m0[0][ii];
				m[1][j] = m0[1][ii];
				m[2][j] = m0[2][ii];
				ii += 1;
			}
				ii = 0;
			for (j = snew - 1; j < ednew; j++){
				m[0][j] = sr3[ii];
				m[1][j] = sc3[ii];
				m[2][j] = sd3[ii];
				ii += 1;
			}
			ed0 = ed; ednew0 = ednew;
		}
		free(rsub0); free(csub0); free(dsub0);
		free(sr3); free(sc3); free(sd3);
	}  
	j = ed0;
	for (i = ednew0; i < len; i++){
		
		m[0][i] = m0[0][j];
		m[1][i] = m0[1][j];
		m[2][i] = m0[2][j];
		j += 1;
	}
	if (inc>0)
	{
		if (dim - inc>0)  
		for (i = 0; i < nsub;i++)
			stnew[i] = stnew[i] + qN;
		else  
		for (i = 0; i < nsub; i++)
			stnew[i] = stnew[i] + hN;
	}
	else 
	for (i = 0; i < nsub; i++)
		stnew[i] = stnew[i] - qN;

	for (i = 0; i < nsub; i++)
		st0[i] = stnew[i];
	for (i = 0; i < 3; i++)                  //释放内存
		free(m0[i]);
	free(m0);
	free(st); free(stnew);
	N0 = len;
	return m;
}