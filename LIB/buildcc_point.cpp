#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE buildcc builds the cartesian coordinates of the tot atoms present
 * in array list according to the reconstruction data.
 ******************************************************************************/
void buildcc_point(FA_Global* FA,float coor[], float dis, float ang, float dih){
  int i;
  //int j;
  float x[4],y[4],z[4];
  float a,b,c,op,cx,cy,cz,d,xn,yn,zn,ct,st,xk,yk,zk;
  
  //printf("%f %f %f\n",dis,ang,dih);
  //PAUSE;


  for(i=1;i<=3;i++){
    x[i]=0.0f+FA->ori[0];
    y[i]=0.0f+FA->ori[1];
    z[i]=0.0f+FA->ori[2];
    if(i==1){x[i]=1.0f+FA->ori[0];}
    if(i==3){y[i]=1.0f+FA->ori[1];}
  }
  
  a=y[1]*(z[2]-z[3])+y[2]*(z[3]-z[1])+y[3]*(z[1]-z[2]);
  b=z[1]*(x[2]-x[3])+z[2]*(x[3]-x[1])+z[3]*(x[1]-x[2]);
  c=x[1]*(y[2]-y[3])+x[2]*(y[3]-y[1])+x[3]*(y[1]-y[2]);
  op=sqrt(a*a+b*b+c*c);
  
  cx=a/op;
  cy=b/op;
  cz=c/op;
  
  a=x[2]-x[1];
  b=y[2]-y[1];
  c=z[2]-z[1];
  
  d=1.0f/sqrt(a*a+b*b+c*c);
  op=dis*d;
  xn=a*op;
  yn=b*op;
  zn=c*op;
  
  a=cx*cx;
  b=cy*cy;
  c=cz*cz;
  
  ct=(float)(cos(ang*PI/180.0f));
  st=(float)(-sin(ang*PI/180.0f));
  op=1.0f-ct;
  
  xk=(cx*cz*op-cy*st)*zn+((1.0f-a)*ct+a)*xn+(cx*cy*op+cz*st)*yn;
  yk=(cy*cx*op-cz*st)*xn+((1.0f-b)*ct+b)*yn+(cy*cz*op+cx*st)*zn;
  zk=(cz*cy*op-cx*st)*yn+((1.0f-c)*ct+c)*zn+(cz*cx*op+cy*st)*xn;
  
  ct=(float)(cos(dih*PI/180.0f));
  st=(float)(sin(dih*PI/180.0f));
  
  cx=(x[2]-x[1])*d;
  cy=(y[2]-y[1])*d;
  cz=(z[2]-z[1])*d;
  
  op=1.0f-ct;
  
  a=cx*cx;
  b=cy*cy;
  c=cz*cz;
  
  x[0]=(cx*cz*op-cy*st)*zk+((1.0f-a)*ct+a)*xk+(cx*cy*op+cz*st)*yk+x[1];
  y[0]=(cy*cx*op-cz*st)*xk+((1.0f-b)*ct+b)*yk+(cy*cz*op+cx*st)*zk+y[1];
  z[0]=(cz*cy*op-cx*st)*yk+((1.0f-c)*ct+c)*zk+(cz*cx*op+cy*st)*xk+z[1];
  
  //printf("%f %f %f\n",x[0],y[0],z[0]);

  *coor=x[0];
  *(coor+1)=y[0];
  *(coor+2)=z[0];
  
  
  return;
}
