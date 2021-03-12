#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"


CR_TEST_RET cr_test(Uint32_Dat* cf,int thd, Uint32_Dat*nc)
{
	CR_TEST_RET  ret;
	SEP sep = separate0(cf->dat, thd, cf->len, 1);
	ret.cf2.dat = sep.rw;
	ret.cf2.len = sep.lrw;

	ret.cf1.dat = cf->dat;
	ret.cf1.len = sep.lrk;

	ret.nc1.dat = (unsigned int *)calloc(thd, sizeof(unsigned int));
	ret.nc1.len = thd;
	for (int i = 0; i < thd; i++)
	{
		ret.nc1.dat[i] = nc->dat[i];
	}

	ret.nc2.len = nc->len - thd;
	ret.nc2.dat = (unsigned int *)calloc(ret.nc2.len, sizeof(unsigned int));
	for (int i = thd; i < nc->len; i++)
	{
		ret.nc2.dat[i - thd] = nc->dat[i];
	}
	ret.z.dat = (unsigned char *)calloc((sep.lensepbit / 8 + 1), sizeof(unsigned char));
	ret.z.len = sep.lensepbit;
	Uint8_Dat z;
	z.dat = sep.sep;
	z.len = sep.lensepbit;

	for (int i = 0; i < (((z.len) >> 3) + 1); i++) //将1变为0，将零变为1.
	{
		ret.z.dat[i] = z.dat[i];
	}

	int cr_count, cri_count;
	float cri;
	k_criterion((void*)(&z), &(ret.cr), &cr_count, &(ret.r), &(ret.lB), &(ret.pr), &(ret.nr), NULL,1);
	//(ret.nr.len)--;//可能是bug
	for (int i = 0; i < ret.nr.len; i++)
	{
		ret.nr.dat[i] = ret.nr.dat[i + 1];
	}
	for (int i = 0; i < (((z.len) >> 3) + 1); i++) //将1变为0，将零变为1.
	{
		(z.dat[i]) ^= 0xff;
	}
	k_criterion((void*)&z, &cri, &cri_count, NULL, NULL, NULL, NULL, NULL, 1);
	ret.sumcr = ret.cr + cri;
	ret.sumcrc = cr_count + cri_count;
	//free(z.dat);
	return ret;
}