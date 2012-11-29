#include "flexaid.h"

int assign_constraint(FA_Global* FA, atom* atoms, resid* residue, constraint* cons)
{
  
  int i,j;
  int k1=0,k2=0;

  int anum,rnum;
  char rnam[4];
  char chn;

  //printf("assigning_constraints.\n");

  rnum=cons->rnum1;
  anum=cons->anum1;
  chn=cons->chn1;
  strcpy(rnam,cons->rnam1);

  //printf("searching %s %c %d %d\n",rnam,chn,rnum,anum);

  for(i=1;i<=FA->res_cnt;i++){
    if(residue[i].number == rnum        &&
       residue[i].chn  == chn           &&
       !strcmp(residue[i].name,rnam)){
      
      for(j=residue[i].fatm[0];j<=residue[i].latm[0];j++){
	if(atoms[j].number == anum) { k1 = j; break; }
      }

      break;
    }
  }

  //printf("found k1: %d\n",k1);

  rnum=cons->rnum2;
  anum=cons->anum2;
  chn=cons->chn2;
  strcpy(rnam,cons->rnam2);

  //printf("searching %s %c %d %d\n",rnam,chn,rnum,anum);

  for(i=1;i<=FA->res_cnt;i++){
    if(residue[i].number == rnum        &&
       residue[i].chn  == chn           &&
       !strcmp(residue[i].name,rnam)){
      
      for(j=residue[i].fatm[0];j<=residue[i].latm[0];j++){
	if(atoms[j].number == anum) { k2 = j; break; }
      }
      
      break;
    }
  }
  
  //printf("found k2: %d\n",k2);

  
  if(k1 != 0 && k2 != 0){
    cons->inum1 = k1;
    cons->inum2 = k2;

    return(1);
  }

  return (0);

}
