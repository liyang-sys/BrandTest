#include <stdio.h>
#include "all.h"
#include "encoding.h"
#include "parameter_setting.h"

int sumBin(Uint8_Dat *binZ)
{
	int sumZ = 0;
	int binBytePos = 0;
	int temp = 0;
	for (int i = 0; i < binZ->len; i++)
	{
		binBytePos = i & 7;
		temp = (((binZ->dat)[i >> 3]) >> (7 - binBytePos)) & 1;
		sumZ += temp;
	}
	return sumZ;
}