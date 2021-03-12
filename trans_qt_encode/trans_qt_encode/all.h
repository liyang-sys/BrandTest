#ifndef ALL_H_INCLUDED
#define ALL_H_INCLUDED

#define uint unsigned int
#define uchar unsigned char

#include <stdlib.h>

/*************************************** r0 ***********************************/

void en_r0 (uint *r, int len_r);

typedef struct {
	uint *ri;
	int lenri;
} KIR;
KIR KW_inversion_runs (uint *r, int lenr);

typedef struct {
	//uchar code;
	uint code;
	int lb;
} SFC;
SFC SFcode (uint sym, uint Nsym);

typedef struct {
	uchar *z;
	int lenzbit;
} GOLINV;
GOLINV GolombInv (uint *r, int lb, int lenr);

typedef struct {
    int nrs;
    int lenres;
} EBNL;

EBNL encode_block_nrs_lenres (uchar *z, int lenzbit, int j);

void encode_block3 (uchar *z, int lenzbit);

void encode_block5 (uchar *z, int lenzbit);

typedef struct {
   int codebook;
   int m, m1, m2, cls, k;
} ESSC;

ESSC encode_stationary_source_cbook(float pw);

void encode_stationary_source_bin(uint *r, int len_r, float pw, int k, int m, int m1, int m2, int cls);

void encode_stationary_source_r_invr (uint *r, int inv, int lastRun, float pw, int lenr, int k, int m, int m1, int m2, int cls);

typedef struct {
	uint *r;
	int lenr;
} DER;
DER de_r0 (uchar *bin, int lenr, int lenbinbit);

typedef struct {
	uint sym;
} DES;
DES deSFcode(uchar *bin, int Nsym);

typedef struct {
	uint *r;
	int lenr;
	int lastBit;
	float p;
} GOL;
GOL Golomb0 (uchar *z, int lenzbit);

typedef struct {
	uchar *z;
	uint *r;
	int lenzbit;
	int lenr;
} DEC;
DEC decode_block3_n0s (uchar *bin, int n0s0, int lenbinbit);

DEC decode_block3_n1s (uchar *bin, int n1s0, int lenbinbit);

DEC decode_block5_n0s (uchar *bin, int n0s0, int lenbinbit);

DEC decode_block5_n1s (uchar *bin, int n1s0, int lenbinbit);

DEC decode_stationary_source_lenr (uchar *bin, int codebook, int lenr0);

DEC decode_stationary_source_leninvr (uchar *bin, int codebook, int leninvr0, int lenbinbit);

DEC decode_stationary_source_r_invr (uchar *bin, int codebook, int lenr, int inv, int lenbinbit);

int sum (uint *r, int lenr);

/*************************************** r1 ***********************************/
//±‡¬Î
int en_r1 (uint *r, float p, int lenr, int *n, int on_off);

int* hist(uint* r,int lenr);

int find_thd(int lenr, int* n, float p);

typedef struct
{
	int thd;
	int idx;
} THE;
THE the_index(int thd);

//typedef struct
//{
//   int lrk;
//   uint *rw;
//   int lrw;
//   uchar *sep;
//   int lensepbit;
//} SEP;
typedef struct
{
	uint *rw;
	uint *rk;
	uint *nk;
	uchar *sep;
	int lrw;
	int lrk;
	int  lnk;
	int lensepbit;
	int thd;
} SEP;
SEP separate0(int *r, int thd, int len, int flg);

int en_runs_sep_sub (int thd, uint *rk, uint *nk, uint *rw, int lenrk, int lenrw);

void en_KsideSub_n1s (uchar *z, int lenzbit, int lenr);

void encode_stationary_source_Nsym (uint *r, int lenr, float p_or_codebook, int Nsym, int on_off, int flg);

void encode_stationary_source_Nsym_sub (uint *r, int len_r, float pw, int Nsym, int on_off, int flg);

void en_KsideSub (uchar *c, int len, int wuns);

void en_Kside(int *rk, int lenrk, int thd);

//“Î¬Î

//uint *de_r1 (uchar *biny, int lenbinbit);
DEC de_r1(uchar *biny, int lenbinbit);

int thd_index_inv(int idx);

DEC decode_stationary_source_lenr (uchar *bin, int codebook, int lenr0);

DEC decode_stationary_source_sumr (uchar *bin, int codebook, int sumr);

DEC decode_block3 (uchar *bin, int lenzbit);

DEC decode_block5 (uchar *bin, int lenzbit);

typedef struct
{
   uint *rk;
   uint *rw;
   int lenrk;
   int lenrw;
} DRSS;
DRSS de_runs_sep_sub (uchar *biny, int thd, int lenrk, int lenrw, int lenbinbit);

DEC de_Kside(uchar *bin, int lenr, int thd, int lenb);

DEC de_KsideSub(uchar *bin, int len);

DEC de_KsideSub_n1s (uchar *bin, int n1s0, int lenbinbit);

DEC decode_stationary_source_Nsym_lenr(uchar *bin, int codebook, int lenr0, int Nsym, int lenbinbit);

uint *separate_inv(uchar *sep, int *rw, int *rk, int len, int thd);

int len_bit_to_byte (int len_bit);

#endif // ALL_H_INCLUDED
