#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE buildic calculates the internal coordinates associated to each
 * atom that is reconstructable on residue gresnum.
 ******************************************************************************/
void buildic_point(FA_Global* FA,float coor[],float *dis, float *ang, float *dih){
  float ca[3],cb[3],cc[3];


  ca[0]=1.0f+FA->ori[0];
  ca[1]=0.0f+FA->ori[1];
  ca[2]=0.0f+FA->ori[2];
  
  cb[0]=0.0f+FA->ori[0];
  cb[1]=0.0f+FA->ori[1];
  cb[2]=0.0f+FA->ori[2];
  
  cc[0]=0.0f+FA->ori[0];
  cc[1]=1.0f+FA->ori[1];
  cc[2]=0.0f+FA->ori[2];
  
  
  //printf("i=%d a=%d b=%d c=%d\n",i,a,b,c);
  
  *dis=distance(coor,ca);
  *ang=bndang(coor,ca,cb);
  *dih=dihedral(coor,ca,cb,cc);
  //printf("%f %f %f\n",*dis,*ang,*dih);
  //PAUSE;

  return;
}
