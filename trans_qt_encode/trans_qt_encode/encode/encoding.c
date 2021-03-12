//#include "encoding.h"
//
// /**
//  * @brief
//  * @param  z[IN]
//  * @param  cr[OUT]
//  * @param  cr_count[OUT] : 若不需要此变量，传递NULL
//  * @param  r[OUT]        : 若不需要此变量，传递NULL
//  * @param  _1B[OUT]      : 若不需要此变量，传递NULL
//  * @param  p[OUT]        : 若不需要此变量，传递NULL
//  * @param  nc[OUT]       : 若不需要此变量，传递NULL
//  * @param  k[OUT]        : 若不需要此变量，传递NULL
//  * @param  flg           : 为0则z=z，为1则1-z （比特流）
//  * @retval none
//  * @note
//  */
//void K_Criterion(Uint8_Dat* z, float* cr, int* cr_count, Uint32_Dat* r, int* lB, float* p, Uint32_Dat* nc, int* k, int flg)
//{
//    GOL  temp;
//    unsigned int N0;
//    Uint8_Dat bak;
//    float q = 0;
//    unsigned int max = 0;
//    int dn = 0;
//    float an = 0;
//
//    int k_flag  = 0;                 // 若k=NULL,则需要在函数里为k动态分配内存，并将k_flag置为1，在退出函数前回收内存
//    int nc_flag = 0;                 // 同上
//    int p_flag  = 0;                 // 同上
//    int r_flag  = 0;
//    int cr_count_flag = 0;
//
//    /* 如果不需要这些参数的返回值，则传参传NULL，但是计算仍需要用到，所以临时为其分配内存 */
//    if(k == NULL){
//       k_flag = 1;
//       k = (int *)malloc(sizeof(int));
//       if(k == NULL){
//          STDOUT_ERROR("k内存分配失败，程序退出.\n");
//          exit(-1);
//       }
//    }
//    if(nc == NULL){
//        nc_flag = 1;
//        nc = (Uint32_Dat *)malloc(sizeof(Uint32_Dat));
//        if(nc == NULL){
//           STDOUT_ERROR("nc内存分配失败，程序退出.\n");
//           exit(-1);
//        }
//    }
//    if(p == NULL){
//        p_flag = 1;
//        p = (float *)malloc(sizeof(float));
//        if(p == NULL){
//           STDOUT_ERROR("p内存分配失败，程序退出.\n");
//           exit(-1);
//        }
//    }
//    if(r == NULL){
//        r_flag = 1;
//        r = (Uint32_Dat*)malloc(sizeof(Uint32_Dat));
//    }
//    if(cr_count == NULL){
//       cr_count_flag = 1;
//       cr_count = (int *)malloc(sizeof(int));
//    }
//
//    /* 实现1-z的操作，并且保存原来的z，在函数退出前恢复原数据 */
//    if(flg){
//       bak.len = z->len;
//       bak.dat = (unsigned char*)malloc(sizeof(unsigned char)*bak.len);
//       memcpy(bak.dat, z->dat, sizeof(unsigned char)*(z->len/8+1) );
//       for(int i=0; i<z->len/8+1; i++){
//         z->dat[i] ^= 0xFF;
//       }
//    }
//
//#if 0
//MATLAB代码：
//if max(z)<=1
//   [r,lB,p]=Golomb0(z); N0=length(r);
//else
//   r=z; N0=length(r); p=N0/sum(r); lB=0;
//end
//
//有点迷，看不懂，z1传进来的是二进制码流，就0和1，max(z)一定会小于等于1的，不知道else有什么用，
//故只实现 [r,lB,p]=Golomb0(z); N0=length(r);
//#endif
//
//    /* [r,lB,p]=Golomb0(z); */
//    temp     =  Golomb0 (z->dat, z->len);
//    r -> dat = temp.r;
//    r -> len = temp.lenr;
//    *lB      = temp.lastBit;
//    *p       = temp.p;
//    /* N0=length(r); */
//    N0       = r->len;
//
//   if(flg){
//      memcpy(z->dat, bak.dat, sizeof(unsigned char)*(z->len/8+1));
//      free(bak.dat);
//      bak.dat = NULL;
//   }
//
///* 调试代码段 */
//#if 0
//STDOUT_DEBUG("lB:%d\n", *lB);
//STDOUT_DEBUG("p:%f\n", *p);
//STDOUT_DEBUG("N0:%d\n", N0);
//#endif
//
//    q = 1 - *p;
//
//    /* 一些辅助变量，因为log函数需要的时间比较久，所以有些固定参数的值就求出来存放在一遍就行了 */
//    /* 减少调用log函数的次数 */
//    float log_q = log(q);
//
//    /* max(r) */
//    max = 0;
//    for(int i=0; i<r->len; i++){
//        if(max < r->dat[i]) max = r->dat[i];
//    }
//
//    /* nc=hist(r, [1:max(r)]); */
//    nc->len = max;
//    nc->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc->len);
//    memset(nc->dat, 0, sizeof(unsigned int)*nc->len);
//    for(int i=0; i<r->len; i++){
//        nc->dat[(unsigned int)(r->dat[i]-1)] += 1;
//    }
//
// /* 调试代码段 */
//#if 0
//int sum = 0;
//STDOUT_DEBUG("nc数组输出\n");
//for(int i=0; i<nc->len; i++){
//    STDOUT_DEBUG("%4d:%2d\n", i+1,nc->dat[i]);
//    sum += nc->dat[i];
//}
//STDOUT_DEBUG("长度：%d 和为%d\n", nc->len, sum);
//while(1);
//#endif
//
//    /* k=floor(-log(N0)/log(q)); */
//    *k = floor( -log(N0)/log_q );
//
//    /* cr_count=sum(nc(k+1:length(nc))); */
//    *cr_count = 0;
//    for(int i=*k; i<nc->len; i++)
//        *cr_count += nc->dat[i];
//
//    *cr = 0;
//    dn = 0;
//    for(int n=*k; n<nc->len; n++){
//        if(nc->dat[n] > 0){
//           if(dn == 0){
//              an = (float)nc->dat[n]/N0;
//              *cr += log(an/(*p)) - n*log_q;
//              //STDOUT_DEBUG("(1)cr = %f\n", *cr);
//           }else{
//              an = ((float)nc->dat[n]/N0) / (dn+1);
//              *cr += log(an/(*p)) - (n - ceil((float)dn/2) )*log_q;
//              //STDOUT_DEBUG("(2)cr = %f dn=%d log(q)=%f\n", *cr, dn, log_q);
//           }
//           dn = 0;
//        }else{
//           dn += 1;
//        }
//    }
//
// /* 调试代码段 */
//#if 0
//STDOUT_DEBUG("cr = %f\n", *cr);
//STDOUT_DEBUG("cr_count = %d\n", *cr_count);
//#endif
//
//    /*! 根据标志位回收内存 */
//    if(k_flag){
//        free(k);
//    }
//    if(nc_flag){
//        free(nc->dat);
//        free(nc);
//    }
//    if(p_flag){
//        free(p);
//    }
//    if(r_flag){
//        free(r->dat);
//        free(r);
//    }
//    if(cr_count_flag){
//        free(cr_count);
//    }
//
//    return;
//}
//
// /**
//  * @brief
//  * @param  n[IN]
//  * @param  p[IN]
//  * @param  thd_arr[IN]
//  * @param  idx[OUT]
//  * @retval none
//  * @note
//  */
//void Find_The_Idx_Round(Uint32_Dat* n, float p, int * thd_arr, int *idx)
//{
//   unsigned int sm   = 0;
//   float        psm  = 0;
//   int          lenn = 0;
//   unsigned int ssm  = 0;
//   unsigned int sum = 0, sum2 = 0;
//
//   for(int i=0; i<n->len; i++ ){
//      sm += n->dat[i];
//      lenn ++;
//   }
//   psm  = p*sm;
//   *idx  = 1;
//   ssm  = n->dat[0];
//   lenn = lenn;
//
///* 调试代码段 */
//#if 0
//STDOUT_INFO("psm  = %f\n", psm);
//STDOUT_INFO("ssm  = %d\n", ssm);
//STDOUT_INFO("lenn = %d\n", lenn);
//#endif
//   /* Hint: thd_arr C语言和MATLAB并没有实现索引同步，C语言指针索引减1 */
//   while( ssm<psm && (*idx)<24 && thd_arr[(*idx)]<lenn ){
//      /* sum(n(thd_arr(idx)+1:thd_arr(idx+1))); */
//      sum = 0;
//      for(int i=thd_arr[(*idx)-1]+1; i<thd_arr[(*idx)]; i++){
//        sum += n->dat[i-1];
//      }
//      ssm += sum;
//      (*idx) ++;
//   }
//
///* 调试代码段 */
//#if 0
//STDOUT_INFO("ssm=%d\n", ssm);
//STDOUT_INFO("idx=%d\n", (*idx));
//#endif
//
//   if((*idx)>1){
//      /* sum(n(1:thd_arr(idx-1))) */
//      sum = 0;
//      for(int i=0; i<thd_arr[(*idx)-2]; i++)
//        sum += n->dat[i-1];
//      /* sum( n(1:thd_arr(idx))) */
//      sum2 = 0;
//      for(int i=0; i<thd_arr[(*idx)-1]; i++)
//        sum2 += n->dat[i-1];
//      /* if psm-sum(n(1:thd_arr(idx-1)))  <  sum( n(1:thd_arr(idx)))-psm */
//      if(psm-sum < sum2-psm){
//        (*idx) --;
//      }
//   }
//   //STDOUT_INFO("idx=%d\n", *idx);
//
//   return;
//}
//
// /**
//  * @brief
//  * @param  r[IN]
//  * @param  thd[IN]
//  * @param  sep[OUT]
//  * @param  rw[OUT]
//  * @param  rk[OUT]
//  * @retval none
//  * @note
//  */
//void Separate(Int32_Dat *cf0, int thd, Uint8_Dat* sep, Uint32_Dat *rw, Uint32_Dat *rk )
//{
//   SEP  temp;
//
//   temp = separate0(cf0->dat, thd, cf0->len,  0);
//
//   rw->len = temp.lrw;
//   rw->dat = temp.rw;
//   rk->len = temp.lrk;
//   sep->dat = temp.sep;
//   sep->len = temp.lensepbit;
//
///* 调试代码段 */
//#if 0
//   STDOUT_DEBUG("rw len = %d\n", rw->len);
//   STDOUT_DEBUG("rk len = %d\n", rk->len);
//   STDOUT_DEBUG("\n");
//   STDOUT_DEBUG("sp len = %d\n", sep->len);
//   STDOUT_DEBUG("SEP:\n");
//   for(int i=0; i < sep->len/8 + 1; i++){
//      printf("  %4d:0x%2X", i+1 , sep->dat[i]);
//      if((i+1)%8 == 0) printf("\n");
//   }
//   STDOUT_DEBUG("\n");
//#endif
//   return;
//}
//
// /**
//  * @brief
//  * @param  cf0[IN]
//  * @param  thd_arr[IN]
//  * @param  nc[IN]
//  * @param  idx0[IN]
//  * @param  idx[OUT]
//  * @param  r[OUT]
//  * @param  lB[OUT]
//  * @param  pr[OUT]
//  * @param  nr[OUT]
//  * @param  sumcr[OUT]
//  * @param  sumcrc[OUT]
//  * @param  cf1k[OUT]
//  * @param  nc1[OUT]
//  * @param  nc1k[OUT]
//  * @param  z1[OUT]
//  * @retval none
//  * @note
//  */
//void Find_Thd4subs(Int32_Dat* cf0, Int32_Dat* thd_arr, Uint32_Dat* nc, int idx0, int *idx,  Uint32_Dat* r, int *lB,
//                   float *pr, Uint32_Dat* nr, float* sumcr, float* sumcrc, Uint32_Dat* cf1,  Uint32_Dat* cf1k,
//                   Uint32_Dat* nc1, Uint32_Dat* nc1k, Uint8_Dat* z1)
//{
//    float p = 0.5;
//    *idx = 0;
//    int   thd  = 0;
//    float cr   = 0;
//    int   crc  = 0;
//    float cri  = 0;
//    int   crci = 0;
//    int   lBi  = 0;
//    Uint32_Dat ri;
//    unsigned int maxcr = 0;
//    float sumcr0   = 0;
//    float sumcrc0  = 0;
//    int   flg      = 0;
//    int   max_idx  = 0;
//    int   sum_nc1k = 0;
//    int   sum_nc   = 0;
//
//    Find_The_Idx_Round(nc, p, thd_arr->dat, idx);
//
//    if(idx0>*idx && thd_arr->dat[*idx]<nc->len){
//        (*idx) ++;
//    }
//
//    thd = thd_arr->dat[*idx-1];
//
//    Separate(cf0, thd, z1, cf1, cf1k);
//
//    /* cf1=cf1-thd;  */
//    for(int i=0; i<cf1->len; i++){
//        cf1->dat[i] += -thd;
//    }
//
//
//    /* nc1k=nc(1:thd); */
//    nc1k->len = thd;
//    nc1k->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1k->len);
//    memcpy(nc1k->dat, nc->dat, sizeof(unsigned int)*nc1k->len);
//
//    /* nc1=nc(thd+1:length(nc)); */
//    nc1->len = nc->len - thd;
//    nc1->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1->len);
//    memcpy(nc1->dat, &(nc->dat[thd]), sizeof(unsigned int)*nc1->len);
//
// /* 调试代码段 */
//#if 0
//STDOUT_DEBUG("nc:\n");
//for(int i=0; i<nc->len; i++)
//    printf("  %d", nc->dat[i]);
//printf("\n");
//STDOUT_DEBUG("nc1k:\n");
//for(int i=0; i<nc1k->len; i++)
//    printf("  %d", nc1k->dat[i]);
//printf("\n");
//STDOUT_DEBUG("nc1:\n");
//for(int i=0; i<nc1->len; i++)
//    printf("  %d", nc1->dat[i]);
//printf("\n");
//while(1);
//#endif
//
//    /* [cr,crc,r,lB,pr,nr]=k_criterion(z1); */
//    K_Criterion(z1, &cr, &crc, r, lB, pr, nr, NULL, 0);
//    /* [cri,crci,ri,lBi]=k_criterion(1-z1); */
//    K_Criterion(z1, &cri, &crci, &ri, &lBi, NULL, NULL, NULL, 1);
//
//    /* maxcr=max([cr cri]); */
//    maxcr = cr>cri ? cr : cri;
//    maxcr = maxcr;                // 为了清除编译器unused的警告
//
//    sumcr0 = 0;
//    *sumcr  = cr  + cri;
//    *sumcrc = crc + crci;
//
//    flg = 0;
//    /* max_idx=numel(thd_arr); */
//    max_idx = thd_arr->len;
//
//    /* sum(nc) */
//    sum_nc = 0;
//    for(int i=0; i<nc->len; i++)
//        sum_nc += nc->dat[i];
//    /* sum(nc1k) */
//    sum_nc1k = 0;
//    for(int i=0; i<nc1k->len; i++)
//        sum_nc1k += nc1k->dat[i];
//    while( *sumcr>sumcr0 && (float)sum_nc1k/sum_nc < 0.86 && *idx<max_idx && thd_arr->dat[*idx]<nc->len){
//        flg = 1;
//        sumcr0  = *sumcr;
//        /* 好像这句话没有什么作用，但还是保留下来吧 */
//        sumcrc0 = *sumcrc;
//        sumcrc0 = sumcrc0;
//        *idx += 1;
//        thd   = thd_arr->dat[(unsigned int)(*idx-1)];
//        /* [z1,cf1,cf1k]=separate(cf0,thd); */
//        Separate(cf0, thd, z1, cf1, cf1k);
//       /* cf1=cf1-thd;  */
//       for(int i=0; i<cf1->len; i++){
//         cf1->dat[i] += -thd;
//       }
//
//        /* nc1k=nc(1:thd); */
//        nc1k->len = thd;
//        free(nc1k->dat);
//        nc1k->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1k->len);
//        memcpy(nc1k->dat, nc->dat, sizeof(unsigned int)*nc1k->len);
//        /* nc1=nc(thd+1:length(nc)); */
//        nc1->len = nc->len - thd;
//        free(nc1->dat);
//        nc1->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1->len);
//        memcpy(nc1->dat, &(nc->dat[thd]), sizeof(unsigned int)*nc1->len);
//
//       /* [cr,crc,r,lB,pr,nr]=k_criterion(z1); */
//       K_Criterion(z1, &cr, &crc, r, lB, pr, nr, NULL, 0);
//       /* [cri,crci,ri,lBi]=k_criterion(1-z1); */
//       K_Criterion(z1, &cri, &crci, &ri, &lBi, NULL, NULL, NULL, 1);
//       /* maxcr=max([cr cri]); */
//       maxcr = cr>cri ? cr : cri;
//       *sumcr  = cr  + cri;
//       *sumcrc = crc + crci;
//    }
//
//    if(flg == 1){
//       *idx += -1;
//       thd   = thd_arr->dat[(unsigned int)(*idx-1)];
//       /* [z1,cf1,cf1k]=separate(cf0,thd); */
//        Separate(cf0, thd, z1, cf1, cf1k);
//       /* cf1=cf1-thd;  */
//       for(int i=0; i<cf1->len; i++){
//         cf1->dat[i] += -thd;
//       }
//
//        /* nc1k=nc(1:thd); */
//        nc1k->len = thd;
//        free(nc1k->dat);
//        nc1k->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1k->len);
//        memcpy(nc1k->dat, nc->dat, sizeof(unsigned int)*nc1k->len);
//        /* nc1=nc(thd+1:length(nc)); */
//        nc1->len = nc->len - thd;
//        free(nc1->dat);
//        nc1->dat = (unsigned int *)malloc(sizeof(unsigned int)*nc1->len);
//        memcpy(nc1->dat, &(nc->dat[thd]), sizeof(unsigned int)*nc1->len);
//
//       /* [cr,crc,r,lB,pr,nr]=k_criterion(z1); */
//       K_Criterion(z1, &cr, &crc, r, lB, pr, nr, NULL, 0);
//       /* [cri,crci,ri,lBi]=k_criterion(1-z1); */
//       K_Criterion(z1, &cri, &crci, &ri, &lBi, NULL, NULL, NULL, 1);
//       /* maxcr=max([cr cri]); */
//       maxcr = cr>cri ? cr : cri;
//       *sumcr  = cr  + cri;
//       *sumcrc = crc + crci;
//    }
//
//    free(ri.dat);
//    ri.dat = NULL;
///* 调试代码段 */
//#if 1
//STDOUT_DEBUG("idx = %d\n", *idx);
//STDOUT_DEBUG("lB = %d\n", *lB);
//STDOUT_DEBUG("pr = %f\n", *pr);
//STDOUT_DEBUG("nr len=%d\n", nr->len);
//STDOUT_DEBUG("sumcr = %f\n", *sumcr);
//STDOUT_DEBUG("sumcrc = %f\n", *sumcrc);
//#endif
//
//    return;
//}
//
// /**
//  * @brief
//  * @param  cf0[IN]       : 量化的得到的一维数组
//  * @param  sn[IN]        :
//  * @param  max_cf0[IN]   : cf0的最大值
//  * @retval none
//  * @note
//  */
//void En_Sub3d_Sub2(Int32_Dat *cf0, Uint32_Dat *sn, char maxcf0)
//{
//    unsigned char maxcf = 0;
//    float  p     = 0;
//    float  p1    = 0;
//    float  p2    = 0;
//    int Nidx     = 0;
//    int lenw     = 0;
//    int lenc     = 0;
//    int iter     = 0;
//    int idx0     = 0;
//    int crc      = 20;
//
//    /*! Find_Thd4subs 需要传递回来的参数 */
//    int thdidx   = 0;
//    Uint32_Dat r;
//    int lB = 0;
//    float pr = 0;
//    Uint32_Dat nr;
//    float scr = 0;
//    float scrc = 0;
//    Uint32_Dat  cf1;
//    Uint32_Dat  cf1k;
//    Uint32_Dat nc1;
//    Uint32_Dat nc1k;
//    Uint8_Dat    z;
//
//
//    /* cf0=cf0+1 */
//    for(int i=0; i<cf0->len; i++){
//        cf0->dat[i] += 1;
//    }
//
//    maxcf =  maxcf0 + 1;
//
//    Uint32_Dat nc;
//    nc.len = maxcf;
//    nc.dat = (unsigned int *)malloc(sizeof(unsigned int)*nc.len);
//    if(nc.dat == NULL){
//        printf("nc内存分配失败，程序退出.\n");
//        exit(-1);
//    }
//    memset(nc.dat, 0, sizeof(unsigned int)*nc.len);
//
//    for(int i=0; i<cf0->len; i++){
//        nc.dat[(unsigned char)(cf0->dat[i] - 1)] += 1;
//    }
//
///* 调试代码段 */
//#if 0
//    STDOUT_INFO("nc数组输出:\n");
//    for(int i=0; i<nc.len; i++)
//        printf("  %2d", nc.dat[i]);
//    printf("\n");
//#endif
//
//#if 0
//MATLAB代码段:
//th=[cumprod(2*ones(1,11)); 3*cumprod([1 2*ones(1,10)])]; th=reshape(th,1,22); thd_arr=[1 th];
//Nidx=length(thd_arr);
//C语言手动生成
//#endif
//
//    Int32_Dat thd_arr;
//    thd_arr.len =  23;
//    thd_arr.dat = (int *)malloc(sizeof(int)*thd_arr.len);
//    thd_arr.dat[0]  = 1;
//    thd_arr.dat[1]  = 2;
//    thd_arr.dat[2]  = 3;
//    thd_arr.dat[3]  = 4;
//    thd_arr.dat[4]  = 6;
//    thd_arr.dat[5]  = 8;
//    thd_arr.dat[6]  = 12;
//    thd_arr.dat[7]  = 16;
//    thd_arr.dat[8]  = 24;
//    thd_arr.dat[9]  = 32;
//    thd_arr.dat[10] = 48;
//    thd_arr.dat[11] = 64;
//    thd_arr.dat[12] = 96;
//    thd_arr.dat[13] = 128;
//    thd_arr.dat[14] = 192;
//    thd_arr.dat[15] = 256;
//    thd_arr.dat[16] = 384;
//    thd_arr.dat[17] = 512;
//    thd_arr.dat[18] = 768;
//    thd_arr.dat[19] = 1024;
//    thd_arr.dat[20] = 1536;
//    thd_arr.dat[21] = 2048;
//    thd_arr.dat[22] = 3072;
//    Nidx = 23;
//
//    p  = 0.7;
//    p1 = 0.7;
//    p2 = 0.8;
//
//    lenw = cf0->len;
//    lenc = 1600;
//
//    iter = 1;
//    Int32_Dat  cf;
//    cf.len = cf0->len;
//    cf.dat = (int *)malloc(sizeof(int)*cf0->len);
//    if(cf.dat == NULL){
//        printf("cf内存分配失败，程序退出.\n");
//        exit(-1);
//    }
//
//    /* cf=cf0; */
//    memcpy(cf.dat, cf0->dat, sizeof(int)*cf0->len);
//
//    idx0 = 1;
//    /* maxcf=max(cf);  */
//    maxcf = 0;
//    for(int i=0; i<cf.len; i++){
//        if(cf.dat[i]>maxcf) maxcf = cf.dat[i];
//    }
//
//    /* while (lenw>lenc || crc>=12)  && iter<=8  &&  maxcf>=2 */
//    while( (lenw>lenc || crc>=12) && iter<=8 && maxcf>=2 ){
//
///* 调试代码段 */
//#if 0
//STDOUT_INFO("lenw  = %d\n", lenw);
//STDOUT_INFO("lenc  = %d\n", lenc);
//STDOUT_INFO("crc   = %d\n", crc);
//STDOUT_INFO("iter  = %d\n", iter);
//STDOUT_INFO("maxcf = %d\n", maxcf);
//while(1);
//#endif
//      Find_Thd4subs(&cf, &thd_arr, &nc, idx0, &thdidx, &r, &lB, &pr, &nr, &scr, &scrc, &cf1, &cf1k, &nc1 ,&nc1k, &z );
//      STDOUT_INFO("写到这里：\n File:%s\n Line:%d\n", __FILE__, __LINE__);
//      exit(-1);
//    }
//
//    free(nc.dat);
//    nc.dat = NULL;
//    free(cf.dat);
//    cf.dat = NULL;
//    free(cf1.dat);
//    cf1.dat = NULL;
//    free(cf1k.dat);
//    cf1k.dat = NULL;
//
//    return;
//}
