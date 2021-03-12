#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char




/** @brief
  * @param  r[IN]             : 会被修改
  * @param  kcr[IN]
  * @param  sep[OUT]
  * @retval none
  * @note   对应关系表:
  *         1 - sep{i,1}      (*sep)[i].sep  (*sep)[i].lensepbit
  *         2 - sep{i,2}      (*sep)[i].rk   (*sep)[i].lrk
  *         3 - sep{i,3}      (*sep)[i].thd
  *         4 - sep{i,4}      (*sep)[i].nk   (*sep)[i].lrk
  */
void en_zone_sub_new(Uint32_Dat* r, float kcr, SEP **sep)
{
	union data {
		uint a;
		uchar b[4];
	} rem;
	extern int ptr;
	extern uchar *bin;
	int x;
	SFC sfc1;
	GOL gol;

	float p1 = 0, p2 = 0;
	float p[7] = { 0 };
	int lenr = 0;
	int nsep = 0;
	int* n = NULL;
	int * nw = NULL;
	int * rw = NULL;
	int * nk = NULL;
	int lenw = 0;
	int idx = 0;
	int thd = 0;
	unsigned char cbook = 0;
	float prob = 0;


	/* p1=1-1/4-1/16+1/64; */
	p1 = 0.7031;
	/* p2=1-1/4+1/16; */
	p2 = 0.8125;
	/* lenr=length(r); */
	lenr = r->len;
	/* if kcr<=60  |  (kcr>60 & lenr<1000) */
	if (kcr <= 60 || (kcr > 60 && lenr < 1000)) {
		nsep = 2;
		/* p=[p1 p2]; */
		p[1] = p1;
		p[2] = p2;
		/* biny=[biny SFcode(1,5)]; */
		sfc1 = SFcode(1, 5);
		x = ptr & 7;
		rem.a = sfc1.code;
		rem.a = rem.a << (16 - x - sfc1.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc1.lb; x += sfc1.lb; x &= 7;
	}
	else if (kcr < 450) { /* elseif kcr<450 */
		nsep = 3;
		/* p=[p1 p1 p2]; */
		p[1] = p1;
		p[2] = p1;
		p[3] = p2;
		/* biny=[biny SFcode(2,5)]; */
		sfc1 = SFcode(2, 5);
		x = ptr & 7;
		rem.a = sfc1.code;
		rem.a = rem.a << (16 - x - sfc1.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc1.lb; x += sfc1.lb; x &= 7;
	}
	else if (kcr < 3600) { /* elseif kcr<3600 */
		nsep = 4;
		/* p=[p1 p1 p1 p2]; */
		p[1] = p1;
		p[2] = p1;
		p[3] = p1;
		p[4] = p2;
		/* biny=[biny SFcode(3,5)]; */
		sfc1 = SFcode(3, 5);
		x = ptr & 7;
		rem.a = sfc1.code;
		rem.a = rem.a << (16 - x - sfc1.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc1.lb; x += sfc1.lb; x &= 7;
	}
	else if (kcr < 16000) { /* elseif kcr<16000 */
		nsep = 5;
		/* p=[p1 p1 p1 p1 p2]; */
		p[1] = p1;
		p[2] = p1;
		p[3] = p1;
		p[4] = p1;
		p[5] = p2;
		/* biny=[biny SFcode(4,5)]; */
		sfc1 = SFcode(4, 5);
		x = ptr & 7;
		rem.a = sfc1.code;
		rem.a = rem.a << (16 - x - sfc1.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc1.lb; x += sfc1.lb; x &= 7;
	}
	else {
		nsep = 6;
		/* p=[p1 p1 p1 p1 p1 p2]; */
		p[1] = p1;
		p[2] = p1;
		p[3] = p1;
		p[4] = p1;
		p[5] = p1;
		p[6] = p2;
		/* biny=[biny SFcode(5,5)]; */
		sfc1 = SFcode(5, 5);
		x = ptr & 7;
		rem.a = sfc1.code;
		rem.a = rem.a << (16 - x - sfc1.lb);
		bin[ptr >> 3] |= rem.b[1];
		bin[(ptr >> 3) + 1] |= rem.b[0];
		ptr += sfc1.lb; x += sfc1.lb; x &= 7;
	}
	/* nnsep=nsep; 冗余 */

	//printf("nsep = %d\n", nsep);
	//printf("(1)ptr = %d\n", ptr);

	/* %----------------------------------------------------------------------- */
	/* th=[cumprod(2*ones(1,11)); 3*cumprod([1 2*ones(1,10)])]; th=reshape(th,1,22); th=[1 th 4096]; */
	unsigned int th[25] = { 0,1,2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,
							768,1024,1536,2048,3072,4096 };
	/* n=hist(r,[1:max(r)]); */
	n = hist(r->dat, r->len);

	/* rw=r; */
	rw = (int *)r->dat;
	rw = rw;              // 保留，格式统一，实际无意义s
	/* nw=n; */
	nw = n;
	/* lenw=length(nw); */
	int max = 0;
	for (int i = 0; i < r->len; i++)
	{
		if (max < r->dat[i])
		{
			max = r->dat[i];
		}
	}
	lenw = max;
	/* sep=cell(nsep,4); */
	*sep = malloc(sizeof(SEP)*nsep);
	if (sep == NULL) {
		perror("sep");
	}
	/* 调试代码段 */
#if 0
	printf("n len:%d\n", lenw);
	for (int i = 0; i < lenw; i++) {
		printf("(%3d):%2d ", i, n[i + 1]);
		if ((i + 1) % 8 == 0) printf("\n");
	}
#endif

	/* for nsp=1:nsep */
	for (int i = 0; i < nsep; i++) {
		SEP temp;
		/* idx=find_thd_idx_round(nw,p(nsp),th); */
		Int32_Dat temp_nw;
		temp_nw.dat = nw + 1;
		temp_nw.len = lenw;
		find_the_idx_round((Uint32_Dat *)&temp_nw, p[i + 1], (int *)th + 1, &idx);
		//printf("idx = %d\n", idx);
		/* thd=th(idx); */
		thd = th[idx];
		/* bin=encode_stationary_source(idx, 3); biny=[uint8(biny) bin]; */
		encode_stationary_source_bin((unsigned int *)&idx, 1, 3, 0, 0, 0, 0, 0);
		//printf("(1.1)ptr = %d\n", ptr);
		/* [sept,rw,rk]=separate(rw,thd);   rw=rw-thd; */
		temp = separate0((int *)r->dat, thd, r->len, 1);
		/* 调试代码段 */
#if 0
		printf("lenbit sep = %d\n", temp.lensepbit);
		printf("len rk = %d\n", temp.lrk);
		printf("len rw = %d\n", temp.lrw);
#endif
		/* sep{nsp,1}=sept; sep{nsp,2}=rk; */
		temp.rk = r->dat;
		r->dat = temp.rw;
		r->len = temp.lrw;
		temp.rw = NULL;
		temp.lrw = 0;
		/* nk=nw(1:thd); sep{nsp,4}=nk; */
		nk = malloc(sizeof(int)*(thd + 1));
		if (nk == NULL) {
			perror("nk");
			exit(-1);
		}
		memcpy(nk, nw, sizeof(int)*(thd + 1));
		temp.nk = (unsigned int *)nk;
		temp.lnk = thd;
		nk = NULL;
		/* nw=nw(thd+1:lenw); */
		memcpy(nw + 1, nw + thd + 1, sizeof(int)*(lenw - thd));
		/* lenw=lenw-thd; */
		lenw = lenw - thd;
		/* sep{nsp,3}=thd; */
		temp.thd = thd;
		(*sep)[i] = temp;
	}
	//printf("(2)ptr = %d\n", ptr);

/* 调试代码段 */
#if 0
	for (int i = 0; i < nsep; i++) {
		printf("i = %d\n", i + 1);
		printf("lenbit sep = %d\n", (*sep)[i].lensepbit);
		printf("len rk = %d\n", (*sep)[i].lrk);
		printf("thd = %d\n", (*sep)[i].thd);
		printf("len nk = %d\n", (*sep)[i].lnk);
	}
	while (1);
#endif

	/*! %encode sep{nsep,1} ------------------------------------------- */
	/* cbook=17; */
	cbook = 17;
	/*  biny=en_runs_0seps(sep{nsep,1},cbook,biny); */
	Uint8_Dat temp_sep;
	temp_sep.dat = (*sep)[nsep - 1].sep;
	temp_sep.len = (*sep)[nsep - 1].lensepbit;
	en_runs_0seps(&temp_sep, cbook);
	//printf("(3)ptr = %d\n", ptr);

	/*! %encode sep{k,1} ( k from 2 to nsep-1 ): */
	/* for nsp = (nsep-1): -1 : 2 */
	for (int nsp = nsep - 1; nsp >= 2; nsp--) {
		/* rk=sep{nsp,2}; thd=sep{nsp,3}; nk=sep{nsp,4}; */
		/* Hint:存放在*sep里 */
		/* [rs,lastBit,ps]=Golomb0(sep{nsp,1}); */
		gol = Golomb0((*sep)[nsp - 1].sep, (*sep)[nsp - 1].lensepbit);
		/* biny=[biny lastBit]; */
		x = ptr & 7;
		bin[ptr >> 3] |= gol.lastBit << (7 - x);
		ptr++;
		/* cbook=4; */
		cbook = 4;
		/*  biny=en_runs_1sep(rs,p2,cbook,length(sep{nsp,1}),biny);    */
		Uint32_Dat rs;
		rs.dat = gol.r;
		rs.len = gol.lenr;
		en_runs_1sep(&rs, p2, cbook, (*sep)[nsp - 1].lensepbit, NULL);
#if 0
		内存释放：
			1 - gol.r
#endif
			free(gol.r);
	}

	//printf("(3.2)ptr = %d\n", ptr);while(1);


	/*! %encode sep{1,1} ----------------------------------------------- */
	/* biny0=biny; 发现此句无效 */
	/* [rs,lastBit,ps]=Golomb0(sep{1,1}); */
	gol = Golomb0((*sep)[0].sep, (*sep)[0].lensepbit);
	/* biny=[biny lastBit]; */
	x = ptr & 7;
	bin[ptr >> 3] |= gol.lastBit << (7 - x);
	ptr++;
	/* cbook=17; */
	cbook = 17;
	/* biny=en_runs_2seps(rs,p2,length(sep{1,1}),cbook,biny); */
	Uint32_Dat rs;
	rs.dat = gol.r;
	rs.len = gol.lenr;
	en_runs_2seps(&rs, p2, (*sep)[0].lensepbit, cbook, NULL);
	free(gol.r);
	rs.dat = NULL;
	rs.len = 0;
	//printf("(4)ptr = %d\n", ptr);

	/*! %encode rw: ---------------------------------------------------- */
	/* [bin]=en_r0(rw); biny=[biny bin]; */
	en_r0(r->dat, r->len);
	//printf("(4)ptr = %d\n", ptr);while(1);

	/*! %encode sep{k,2}: --------------------------------------------------------- */
   /* for nsp=nsep:-1:1 */
	for (int nsp = nsep; nsp >= 1; nsp--) {
		/* rk=sep{nsp,2}; thd=sep{nsp,3}; nk=sep{nsp,4}; */
		/* Hint:存放在*sep里 */
		thd = (*sep)[nsp - 1].thd;
		Uint32_Dat rk_temp;
		Uint32_Dat nk_temp;
		rk_temp.dat = (unsigned int *)(*sep)[nsp - 1].rk;
		rk_temp.len = (*sep)[nsp - 1].lrk;
		nk_temp.dat = (unsigned int *)(*sep)[nsp - 1].nk + 1;
		nk_temp.len = (*sep)[nsp - 1].lnk;
		if (thd > 1) {
			probility_estimation(&rk_temp, &nk_temp, &prob);
			//printf("prob = %f len rk = %d len nk = %d\n", prob, rk_temp.len, nk_temp.len);
			if (prob > 0.15 && thd <= 9) {
				/* biny = [biny  0 ]; */
				ptr++;
				/* bin=en_Kside_new(rk,thd,nk); biny = [biny  bin]; */
				en_Kside_new((int *)rk_temp.dat, rk_temp.len, (int *)nk_temp.dat, nk_temp.len, thd);
			}
			else {
				//printf("(5)ptr = %d\n", ptr);
				/* biny=[biny 1]; */
				x = ptr & 7;
				bin[ptr >> 3] |= 1 << (7 - x);
				ptr++;
#if 0
				/* biny=[biny codebook]; */
				essc1 = encode_stationary_source_cbook(prob);
				if (thd <= 5) essc1.codebook = 0;
				sfc1 = SFcode(essc1.codebook + 1, 20);
				x = ptr & 7;
				rem.a = sfc1.code;
				rem.a = rem.a << (16 - x - sfc1.lb);
				bin[ptr >> 3] |= rem.b[1];
				bin[(ptr >> 3) + 1] |= rem.b[0];
				ptr += sfc1.lb; x += sfc1.lb; x &= 7;
#endif
				/* 调试代码段 */
#if 0
				printf("**********************************\n");
				printf("prob = %f\n", prob);
				printf("thd = %d\n", thd);
				for (int i = 0; i < rk_temp.len; i++) {
					printf("(%3d):%4d ", i + 1, rk_temp.dat[i]);
					if ((i + 1) % 8 == 0) printf("\n");
				}
				printf("\n");
#endif
				//printf("(6.1)ptr = %d\n", ptr);
				/* [bin,codebook]=encode_stationary_source_Nsym(rk,prob,thd); biny = [biny bin] */
				encode_stationary_source_Nsym(rk_temp.dat, rk_temp.len, prob, thd, 1, 20);
				//printf("(6.2)ptr = %d\n", ptr);
			}
		}
	}
	//printf("(7)ptr = %d\n", ptr);

	return;
}
