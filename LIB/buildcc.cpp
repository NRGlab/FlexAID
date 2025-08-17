#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE buildcc builds the cartesian coordinates of the tot atoms present
 * in array list according to the reconstruction data.
 ******************************************************************************/
void buildcc(FA_Global* FA,atom* atoms,int tot,int list[]){
  int an,i,j;
  float x[4],y[4],z[4];
  float a,b,c,op,cx,cy,cz,d,xn,yn,zn,ct,st,xk,yk,zk,angPI,dihPI;

  for(an=0;an<tot;an++){
    
    for(i=1;i<=3;i++)
    {
      j=atoms[list[an]].rec[i-1];
      //printf("recj[%d]=%d listan[%d]=%d\n",j,atoms[j].number,list[an],atoms[list[an]].number);
      //PAUSE;
      if(j != 0)
      {
      	//x[1],y[1],z[1] are coordinates from rec[0] of atom list[an]
      	//x[2],y[2],z[2] are coordinates from rec[1] of atom list[an]
      	//x[3],y[3],z[3] are coordinates from rec[2] of atom list[an]
      	x[i]=atoms[j].coor[0];
      	y[i]=atoms[j].coor[1];
      	z[i]=atoms[j].coor[2];
      }
      else if(i==1)
      {
      	x[i]=1.0f+FA->ori[0];
      	y[i]=0.0f+FA->ori[1];
      	z[i]=0.0f+FA->ori[2];
      }
      else if(i==3)
      {
      	x[i]=0.0f+FA->ori[0];
      	y[i]=1.0f+FA->ori[1];
      	z[i]=0.0f+FA->ori[2];
      }
      else
      {
      	x[i]=0.0f+FA->ori[0];
      	y[i]=0.0f+FA->ori[1];
      	z[i]=0.0f+FA->ori[2];
      }	

      // perturb atom coordinates
      x[i]+=1e-10f;
      y[i]+=1e-10f;
      z[i]+=1e-10f;
    }
    
    a=y[1]*(z[2]-z[3])+y[2]*(z[3]-z[1])+y[3]*(z[1]-z[2]);
    b=z[1]*(x[2]-x[3])+z[2]*(x[3]-x[1])+z[3]*(x[1]-x[2]);
    c=x[1]*(y[2]-y[3])+x[2]*(y[3]-y[1])+x[3]*(y[1]-y[2]);
    op=sqrt(a*a+b*b+c*c);
    
    /*
      d=x[1]*(y[3]*z[2]-y[2]*z[3])+
      y[1]*(x[2]*z[3]-x[3]*z[2])+
      z[1]*(x[3]*y[2]-x[2]*y[3]);
      printf("d=%f\n",d);
      //PAUSE;
      //if(d <= 0.0){op *= -1.0;}
    */

    cx=a/op;
    cy=b/op;
    cz=c/op;
    //printf("cx=%f cy=%f cz=%f\n",cx,cy,cz);

    a=x[2]-x[1];
    b=y[2]-y[1];
    c=z[2]-z[1];

    d=1.0f/sqrt(a*a+b*b+c*c);
    op=atoms[list[an]].dis*d;
    xn=a*op;
    yn=b*op;
    zn=c*op;
    //printf("d=%f op=%f xn=%f yn=%f zn=%f\n",d,op,xn,yn,zn);

    a=cx*cx;
    b=cy*cy;
    c=cz*cz;

    //printf("ang=%f\n",atoms[list[an]].ang);
    angPI = (float)(atoms[list[an]].ang*PI/180.0f);
    ct=cos(angPI);
    st=-sin(angPI);

    op=1.0f-ct;

    //ct=cos(atoms[list[an]].ang*PI/180.0);
    //st=-sin(atoms[list[an]].ang*PI/180.0);
    
    //printf("ct=%f st=%f op=%f\n",ct,st,op);

    xk=(cx*cz*op-cy*st)*zn+((1.0f-a)*ct+a)*xn+(cx*cy*op+cz*st)*yn;
    yk=(cy*cx*op-cz*st)*xn+((1.0f-b)*ct+b)*yn+(cy*cz*op+cx*st)*zn;
    zk=(cz*cy*op-cx*st)*yn+((1.0f-c)*ct+c)*zn+(cz*cx*op+cy*st)*xn;

    //printf("xk=%f yk=%f zk=%f\n",xk,yk,zk);
    //printf("dih=%f\n",atoms[list[an]].dih);
    
    dihPI = (float)(atoms[list[an]].dih*PI/180.0f);
    ct=cos(dihPI);
    st=sin(dihPI);

    op=1.0f-ct;

    //ct=cos(atoms[list[an]].dih*PI/180.0);
    //st=sin(atoms[list[an]].dih*PI/180.0);

    cx=(x[2]-x[1])*d;
    cy=(y[2]-y[1])*d;
    cz=(z[2]-z[1])*d;   
    
    a=cx*cx;
    b=cy*cy;
    c=cz*cz;
    //printf("a=%f b=%f c=%f\n",a,b,c);

    x[0]=(cx*cz*op-cy*st)*zk+((1.0f-a)*ct+a)*xk+(cx*cy*op+cz*st)*yk+x[1];
    y[0]=(cy*cx*op-cz*st)*xk+((1.0f-b)*ct+b)*yk+(cy*cz*op+cx*st)*zk+y[1];
    z[0]=(cz*cy*op-cx*st)*yk+((1.0f-c)*ct+c)*zk+(cz*cx*op+cy*st)*xk+z[1];
    
    atoms[list[an]].coor[0]=x[0];
    atoms[list[an]].coor[1]=y[0];
    atoms[list[an]].coor[2]=z[0];
   }


  return;
}
