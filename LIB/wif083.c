#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE wif083 - creates an sphere of radius 1.0 centered at 0.0 with
 * point almost equally spaced in the surface, used to calculate the 
 * the surface in contact. 
 *                                       Translated from FORTRAN and simplified
 ******************************************************************************/
void wif083(FA_Global* FA){

	int iqa,iqb,iqc,jj,k;
	//int i,kmin,j,min;
	double xb,yb,ph,ct,st;
	//float mindst,dst,tmp;

	iqa=0;
	iqb=1;

	while(iqb <= MAX_SPHERE_POINTS){
		iqc=iqa+iqb;
		iqa=iqb;
		iqb=iqc;
	}

	if(iqb > MAX_SPHERE_POINTS){
		iqb=iqa;
		iqa=iqc-iqb;
	}

	FA->tspoints=iqb;
	//printf("using %d (not %d) points in sphere\n",iqb,MAX_SPHERE_POINTS);


	jj=0;
	for(k=1;k<=iqb;k++){
		xb=k/(double)iqb;
		jj += iqa;
		if (jj > iqb) jj -= iqb;
		yb=jj/(double)iqb;
		ph=2.0*PI*yb;
		ct=1.0-2.0*xb;
		st=sqrt(1.0-ct*ct);
		FA->sphere[k-1][0]=(float)(st*cos(ph));
		FA->sphere[k-1][1]=(float)(st*sin(ph));
		FA->sphere[k-1][2]=(float)ct;
		//printf("%3d %8.3f %8.3f %8.3f\n",k,FA->sphere[k-1][0],FA->sphere[k-1][1],FA->sphere[k-1][2]);
	}
	//exit(8);
	return;
}
