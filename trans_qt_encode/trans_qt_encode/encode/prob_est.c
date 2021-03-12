//#include <stdio.h>
//#include "all.h"
//#include "parameter_setting.h"
//#include "encoding.h"
//
//
//float prob_est(Int32_Dat*n)
//{
//	float p;
//	int nt ;
//	if ((n->len) < 3)
//	{
//		p = 1 - ((float)(n->dat[1]) / (float)(n->dat[0]));
//
//	}
//	else
//	{
//		if ((n->len) >= 16)
//		{
//			nt = 4;
//		}
//		else
//		{
//			nt = 2;
//		}
//		int a = 0;
//		for (int i = 0; i < nt; i++)
//		{
//			a += (n->dat[i]);
//		}
//		int b = 0;
//		for (int i = nt; i < 2 * nt ; i++)
//		{
//			b += (n->dat[i]);
//		}
//		if (b > a)
//		{
//			int c = b;
//			b = a;
//			a = c;
//		}
//		float dx = (float)(a - b) / (float)a;
//		p = dx / nt + (1 / (2.0*nt))*(1 - 1 / (float)(nt))*(dx *dx);
//	}
//	return p;
//}

#include <stdio.h>
#include "all.h"
#include "encoding.h"

#define uint unsigned int
#define uchar unsigned char

void prob_est(Uint32_Dat* n, float *p)
{
	int nt = 0;
	int a = 0;
	int b = 0;
	int c = 0;
	float dx = 0;

	/* if length(n)<=3 */
	if (n->len <= 3) {
		/* p=1-n(2)/n(1);  */
		*p = 1 - (float)n->dat[1] / n->dat[0];
	}
	else {
		/*! % 2*nt-term estimation (power to the (1/nt) ): */
		if (n->len >= 16) {
			nt = 4;
		}
		else {
			nt = 2;
		}
		/* a=sum(n(1:nt)); */
		a = 0;
		for (int i = 0; i < nt; i++) {
			a += n->dat[i];
		}
		/* b=sum( n(nt+1 : 2*nt) ); */
		b = 0;
		for (int i = nt; i < 2 * nt; i++) {
			b += n->dat[i];
		}
		/* if b>a,   c=b; b=a; a=c;   end */
		if (b > a) {
			c = b;
			b = a;
			a = c;
		}
		/* dx=(a-b)/a; */
		dx = (float)(a - b) / a;
		/* p= dx/nt  +(1/(2*nt))*(1-1/nt)*(dx*dx);  */
		*p = dx / nt + (1 / (float)(2 * nt))*(1 - 1 / (float)nt)*(dx*dx);

	}

	return;
}
