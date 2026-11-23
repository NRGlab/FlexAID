#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE buildic calculates the internal coordinates associated to each
 * atom that is reconstructable on residue gresnum.
 ******************************************************************************/
void buildic(FA_Global* FA,atom* atoms,resid* residue,int rnum){
  int i,a,b,c;
  float ci[3],ca[3],cb[3],cc[3];
  int rot;

  rot=residue[rnum].rot;
  for(i=residue[rnum].fatm[rot];i<=residue[rnum].latm[rot];i++){
    if(atoms[i].recs == 'm'){
      a=atoms[i].rec[0];
      b=atoms[i].rec[1];
      c=atoms[i].rec[2];

      ci[0]=atoms[i].coor[0];
      ci[1]=atoms[i].coor[1];
      ci[2]=atoms[i].coor[2];

      if(a == 0){
	ca[0]=1.0f+FA->ori[0];
	ca[1]=0.0f+FA->ori[1];
	ca[2]=0.0f+FA->ori[2];
      }else{
	ca[0]=atoms[a].coor[0];
	ca[1]=atoms[a].coor[1];
	ca[2]=atoms[a].coor[2];
      }

      if(b == 0){
	cb[0]=0.0f+FA->ori[0];
	cb[1]=0.0f+FA->ori[1];
	cb[2]=0.0f+FA->ori[2];
      }else{
	cb[0]=atoms[b].coor[0];
	cb[1]=atoms[b].coor[1];
	cb[2]=atoms[b].coor[2];
      }

      if(c == 0){
	cc[0]=0.0f+FA->ori[0];
	cc[1]=1.0f+FA->ori[1];
	cc[2]=0.0f+FA->ori[2];
      }else{
	cc[0]=atoms[c].coor[0];
	cc[1]=atoms[c].coor[1];
	cc[2]=atoms[c].coor[2];
      }

      /*
      printf("ci[0]=%f ci[1]=%f ci[2]=%f\n",ci[0],ci[1],ci[2]);
      printf("ca[0]=%f ca[1]=%f ca[2]=%f\n",ca[0],ca[1],ca[2]);
      printf("cb[0]=%f cb[1]=%f cb[2]=%f\n",cb[0],cb[1],cb[2]);
      printf("cc[0]=%f cc[1]=%f cc[2]=%f\n",cc[0],cc[1],cc[2]);
      */
      //printf("i=%d a=%d b=%d c=%d\n",i,a,b,c);

      //printf("From dist=%f\tFrom distance=%f\n",dist(ci,ca),distance(ci,ca));
      atoms[i].dis=distance(ci,ca);
      //printf("From bndang=%f\tFrom angle=%f\n",bndang(ci,ca,cb),angle(ci,ca,cb));
      atoms[i].ang=bndang(ci,ca,cb);
      //printf("From dihang=%f\tFrom dihedral=%f\n",dihang(ci,ca,cb,cc),dihedral(ci,ca,cb,cc));
      atoms[i].dih=dihedral(ci,ca,cb,cc);
      //printf("%d: %f %f %f\n",i,atoms[i].dis,atoms[i].ang,atoms[i].dih);
      //PAUSE;
    }
  }
  return;
}
