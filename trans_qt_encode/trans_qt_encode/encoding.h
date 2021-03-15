/**
  ******************************************************************************
  * @file    encoding.h
  * @author  H.R
  * @version V1.0.0
  * @date    2020-10-18
  * @brief   图像编码相关API
  ******************************************************************************
  * @attention
  *
  ******************************************************************************
  */

#ifndef __ENCODING_H
#define __ENCODING_H

#include "parameter_setting.h"
#include "all.h"
#include <math.h>

void k_criterion(void* z, float* cr, int* cr_count, Uint32_Dat* r, int* lB, float* p, Uint32_Dat* nc, int* k, int sep_flg);
//void Find_The_Idx_Round(Uint32_Dat* n, float p, int * thd_arr, int *idx);
void find_the_idx_round(Uint32_Dat* n, float p, int * thd_arr, int *idx);
//void En_Sub3d_Sub2(Int32_Dat *cf0, Uint32_Dat *sn, char maxcf0);
void Find_The_Idx_Round(Uint32_Dat* n, float p, int * thd_arr, int *idx);
void Separate(Int32_Dat *cf0, int thd, Uint8_Dat* sep, Uint32_Dat *rw, Uint32_Dat *rk );
void Find_Thd4subs(Int32_Dat* cf0, Int32_Dat* thd_arr, Uint32_Dat* nc, int idx0, int *idx,  Uint32_Dat* r, int *lB,
                   float *pr, Uint32_Dat* nr, float* sumcr, float* sumcrc, Uint32_Dat* cf1,  Uint32_Dat* cf1k,
                   Uint32_Dat* nc1, Uint32_Dat* nc1k, Uint8_Dat* z1);
void en_sub3d_sub2_sub(Uint8_Dat* z, Uint32_Dat* r, int lB, float pr, Uint32_Dat* nr, float scr, float scrc, Uint32_Dat* cfk, Uint32_Dat* nck, int thd);
void en_r_cr(Uint32_Dat* r, int lB, Uint32_Dat* nr, float sumcr, float sumcrc, float pr, Uint8_Dat* z);
void en_2n_code_new(Uint32_Dat* c2n, Uint32_Dat* nc,int typ,int maxc);
void en_z0(Uint8_Dat* z,int wuns);//wuns是cw的长度
DEC de_z0(uchar *bin, int len);
void en_2n_code_sub(Uint32_Dat *cww, Uint32_Dat* cwk, int n);//cww或者cwk可以为NULL，当为NULL时，不加cww或者cwk编进去
DEC de_2n_code_new(uchar *bin, int len, int maxc, int typ);
DEC de_2n_code_sub(uchar *bin, int n, Uint8_Dat* z);
DEC de_2n_code_sub1(uchar *bin, int n, int len);

void en_z0_r0(Uint8_Dat* z, Uint32_Dat* runs, int lastBit, int wuns, float p);//wuns是z的和
DEC de_z0_r0(uchar *bin, int len);
void en_zone_sub2(Uint32_Dat*r);//biny得在调用此函数之前编入bin中去

void en_runs_0seps(Uint8_Dat* sep, unsigned char cbook);
void en_runs_1sep(Uint32_Dat* r, float p, unsigned char cbook, int lenz, SEP* sep);
void prob_est(Uint32_Dat* n, float *p);
DEC de_r_cr(uchar *bin, int len, int lenbinbit);
void de_runs_0seps(unsigned char cbook, Uint8_Dat* sep);
void de_runs_1sep(int lenr, unsigned char cbook, Uint32_Dat* r, int lenbinbit);
//DEC de_zone_sub2(uchar *bin, int lastBit, int lenbinbit);//lenbinbit是bin的位数。
DEC de_zone_sub2(uchar *biny, int lastBit, int lenbinbit);
int sumBin(Uint8_Dat *binZ);
void en_zone_new(Uint32_Dat*r, Uint32_Dat*nc);
int  find_thd_idx(Uint32_Dat*n, float p, int *th);
void en_runs_2seps(Uint32_Dat* r, float p, int lenz, unsigned char cbook, SEP* sep);
void en_Kside_new(int *rk, int l_rk, int *nk, int l_nk, int thd);
void probility_estimation(Uint32_Dat* rk, Uint32_Dat* nk, float *prob);
void en_zone_sub_new(Uint32_Dat* r, float kcr, SEP **sep);
void en_runs_2sepsA(Uint32_Dat* r, float p, unsigned char cbook, SEP* sep);
DEC de_zone_new(uchar *bin, int lenbinbit);
void de_runs_2sepsA(unsigned char cbook, Uint32_Dat* r, int lenbinbit);
DEC de_Kside_new(uchar *bin, int lenr, int thd, int lenby);
void de_runs_2seps(int lenr, unsigned char cbook, Uint32_Dat* r, int lenbinbit);
void de_zone_sub_new(uchar *bin, int lenbinbit, Uint32_Dat *rrw, Uint32_Dat *rk);


void en_sub3d_sub2(float* cf0, Uint8_Dat* sn, int lencf0, int maxcf0);

CR_TEST_RET cr_test(Uint32_Dat* cf, int thd, Uint32_Dat*nc);
DE_S_S_SUB de_sub3d_sub2_sub(uchar *biny, int thd, int len, int lenbinbit);
DE_S_SUB de_sub3d_sub2(uchar *bin, int len, int qctr, int lenbinbit);
int bin_add(uchar *bin1, uchar *bin2, int lenbin2bit);
/*void find_thd4subs(Int32_Dat* cf0, Int32_Dat* thd_arr, Uint32_Dat* nc, int idx0, int *idx, Uint32_Dat* r, int *lB,
	float *pr, Uint32_Dat* nr, float* sumcr, float* sumcrc, Uint32_Dat* cf1, Uint32_Dat* cf1k,
	Uint32_Dat* nc1, Uint32_Dat* nc1k, Uint8_Dat* z1);

//测试
void test_en_2n_code_new();
void test_en_z0();
void test_de_z0();
void test_en_2n_code_sub();
void test_de_2n_code_new();
void test_en_and_de_2n_code_new();
void test_de_r_cr();
void writeBin();
void test_en_zone_sub2();
void test_en_sub3d_sub2();
void test_en_sub3d_sub2_sub();
void test_de_sub3d_sub2_sub();
void test_de_sub3d_sub2();
void test_en_zone_sub_new();*/

void writeToDataFile(float* cf0, Uint8_Dat* sn, int lencf0);
void writeBinToFile(unsigned char *bin, int ptr);
void test();
void encode_stationary_source_Nsym_sub_in_en_Kside_new(uint *r, int lenr, float pw, int Nsym, int flg);
void encode_stationary_source_Nsym_in_en_Kside_new(uint *r, int lenr, float p_or_codebook, int Nsym, int flg);

#endif /* __ENCODING_H */
