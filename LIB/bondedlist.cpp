#include "flexaid.h"
#include "boinc.h"

/******************************************************************************
 * SUBROUTINE bondedlist creates the list of atoms nloops bonds away from atom
 * anum, outputs the list in the blist vector. 
 ******************************************************************************/
void bondedlist(atom* atoms,int anum, int nloops, int *nlist_ptr, int* blist, int* nnbr){
  int j,k,l,m,n;
  int lista[MAX_ATM_HET];
  int nlist,nlista,nlist_nbr;
  int lflag;
  int oldn;
  
  
  nlist=0;
  nlist_nbr=0;
  blist[nlist++]=anum;
  nnbr[nlist_nbr++]=0;
  
  oldn=0;
  for(n=1;n<=nloops;n++){
    nlista=0;
    for(j=oldn;j<nlist;j++){
      for(k=1;k<=atoms[blist[j]].bond[0];k++){
	l=atoms[blist[j]].bond[k];
	lflag=0;
	m=0;
	while(lflag==0 && m<nlist){
	  if(blist[m] == l){lflag=1;}
	  m++;
	}
	if(lflag==0){
	  lista[nlista++]=l;
	  nnbr[nlist_nbr++]=n;
	}
      }
    }
    
    oldn=nlist;

    for(j=0;j<nlista;j++){
      blist[nlist++]=lista[j];
    }

  }
  
  *nlist_ptr=nlist;
  
  /*
    printf("atom[%d].bondedlist=",atoms[anum].number);
    for(j=0;j<nlist;j++){
    printf("%5d",atoms[blist[j]].number);
    }
    printf("\n");
    getchar();  
  */
  
  return;
  
}
