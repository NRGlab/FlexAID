#include "flexaid.h"

/////////////////////////////////////////////////////////////////////////////
/////// This function assigns all shift as it depends on the order //////////
/////// of the atoms built. If an atom is part of a shift complex  //////////
/////// with a non-rebuilt atom add atoms to shift                 //////////
/////////////////////////////////////////////////////////////////////////////

void assign_shift(atom* atoms,resid* residue,int rnum, int tot, int *buildlist, int **fdihlist){
  int i = 0;                             // dumb counter
  int j = 0;                             // dumb counter
  int k = 0;                             // dumb counter
  int l = 0;                             // dumb counter
  int flag;                              // flag to check if an atom defines a flex. bond
  int built;                             // flag to determine if an atom was already built
  int rec;                               // atoms' rec[0]
  
  int nshift;                            // number of shift atoms
  int shiftlist[4];                      // list of atoms to shift
  
  
  int* done;                             // list of flags to check if a flexible bond was already done shifted
  
  done = (int*)malloc((residue[rnum].fdih+1)*sizeof(int));

  for(i=1;i<=residue[rnum].fdih;i++)
    done[i]=0;
  
  /*
  for(i=1;i<=residue[rnum].fdih;i++){
    for(j=0;j<3;j++){
      printf("fdihlist[%d][%d]=%d\n",i,j,atoms[fdihlist[i][j]].number);
    }
  }
  */
  
  flag=0;
  built=0;
  //Skip GPA atoms (0,1,2)
  for(i=3;i<tot;i++){
    
    //printf("buildlist[%d]=%6d(%s)\n",i,buildlist[i],atoms[buildlist[i]].name);
    
    //if atom defines a flexible bond
    //check rec[0] of each atom (x)
    //if its bound to a non-built atom or atom x, add shift
    nshift=0;
    shiftlist[nshift++]=buildlist[i];
    
    flag=0;
    for(j=1;j<=residue[rnum].fdih;j++){
      
      if(done[j] == 1){
	//built = 1;
	continue;
      }
      
      for(k=0;k<3;k++){
	
	if(fdihlist[j][k]!=0 && fdihlist[j][k]==buildlist[i]){
	  l=j;
	  flag=1;
	  break;
	}
	
      }
      
      if (flag) { break; }
    }

    // atom defines a flexible bond
    if(flag){
      //printf("atom[%d] defines a flexible bond\n",buildlist[i]);
      
      rec=atoms[buildlist[i]].rec[0];
      for(j=1;j<=atoms[rec].bond[0];j++){
	//skip atom that defines the flexible bond
	if(atoms[rec].bond[j]==buildlist[i]) continue;
	
	/* NEW VERSION */
	for(k=i+1;k<tot;k++)
	  if(buildlist[k]==atoms[rec].bond[j])
	    shiftlist[nshift++]=atoms[rec].bond[j];
	
      }
      
    }
    
    if(nshift > 1) {
      
      done[l] = 1;
      residue[rnum].bond[l] = shiftlist[0];

      /*
      printf("shiftlist=");
      for(k=0;k<nshift;k++) printf("%6d(%s)",atoms[shiftlist[k]].number,atoms[shiftlist[k]].name);
      printf("\n");
      printf("flexible bond %d done\n",l);
      fflush(stdout);
      */

      for(j=0,k=j+1;j<nshift;j++,k++){
	if(k==nshift) k=0;
	atoms[shiftlist[j]].rec[3] = shiftlist[k];
	//printf("atoms[%d].rec[3]=%d\n",atoms[shiftlist[j]].number,atoms[shiftlist[k]].number);
	atoms[shiftlist[k]].shift = atoms[shiftlist[k]].dih - atoms[shiftlist[j]].dih;
	//printf("atoms[%d].shift=%.3f\n",atoms[shiftlist[k]].number,atoms[shiftlist[k]].shift);
	printf("shiftval= %d %d with value %8.3f\n",atoms[shiftlist[j]].number,atoms[shiftlist[k]].number,atoms[shiftlist[k]].shift);
      }
    }
  }
  
  free(done);

  //printf("assign_shift done!\n");
  return;
}
