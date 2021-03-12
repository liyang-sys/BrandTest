#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h" 
/**
  * @brief
  * @param  r[IN]
  * @param  thd[IN]
  * @param  sep[OUT]
  * @param  rw[OUT]
  * @param  rk[OUT]
  * @retval none
  * @note
  */
void Separate(Int32_Dat *cf0, int thd, Uint8_Dat* sep, Uint32_Dat *rw, Uint32_Dat *rk )
{
   SEP  temp;

   temp = separate0(cf0->dat, thd, cf0->len,  0);

   rw->len = temp.lrw;
   rw->dat = temp.rw;
   rk->len = temp.lrk;
   sep->dat = temp.sep;
   sep->len = temp.lensepbit;

/* µ÷ÊÔ´úÂë¶Î */
#if 0
   STDOUT_DEBUG("rw len = %d\n", rw->len);
   STDOUT_DEBUG("rk len = %d\n", rk->len);
   STDOUT_DEBUG("\n");
   STDOUT_DEBUG("sp len = %d\n", sep->len);
   STDOUT_DEBUG("SEP:\n");
   for(int i=0; i < sep->len/8 + 1; i++){
      printf("  %4d:0x%2X", i+1 , sep->dat[i]);
      if((i+1)%8 == 0) printf("\n");
   }
   STDOUT_DEBUG("\n");
#endif
   return;
}