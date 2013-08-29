#include "flexaid.h"
#include "boinc.h"

// the function updates the matrix of bonded atoms
// the matrix starts at [0,0] with residue[ires].fatm
// works/available for rotamers also

void update_bonded(resid* residue, int tot, int nlist, int* list, int* nbr)
{

  int i,j;             // dumb counters
  int fatm,latm;       // first-last atom of residue i

  /**************************************************/
  /*******     allocate memory for matrix ***********/
  /**************************************************/
  
  if(residue->bonded == NULL)
    {
      
      residue->bonded = (int**)malloc(tot*sizeof(int*));
      
      if(residue->bonded == NULL)
	{
	  fprintf(stderr,"ERROR: Could not allocate memory for residue->bonded\n");
	  Terminate(2);
	}
      else
	{
	  for(i=0;i<tot;i++){

	    residue->bonded[i] = (int*)malloc(tot*sizeof(int));

	    if(residue->bonded[i] == NULL)
	      {
		fprintf(stderr,"ERROR: Could not allocate memory for residue->bonded[%d]\n",i);
		Terminate(2);
	      }	    
	    
	  }
	    
      }
     
      // fill matrix with -1 (nothing found up to nloops)
      for(i=0;i<tot;i++)
	for(j=0;j<tot;j++)
	  residue->bonded[i][j] = -1;
	
      /*
	printf("new bonded residue!\n");
	getchar();
      */
    }


  /**************************************************/
  /*******     update content of matrix   ***********/
  /**************************************************/

  /*
    printf("updating...\nmatrix so far...\n");
    for(i=0;i<tot;++i){for(j=0;j<tot;++j){printf("[%d-%d]=%d\n",i,j,residue->bonded[i][j]);}}
  */

  //  printf("will write matrix to %p\n",residue->bonded);

  // assign first atom to start from 0
  fatm = residue->fatm[0];
  latm = residue->latm[0];
  
  //printf("fatm:%d\tlatm:%d\n",fatm,latm);

  // loop through neighbours list
  for(i=0;i<nlist;i++){
    
    /*
    printf("writing %d against %d with value %d\n",list[0]-fatm,list[i]-fatm,nbr[i]);
    printf("internal atom numbers list[0]=%d\tlist[%d]=%d\tfatm=%d\n",
    	   list[0],i,list[i],fatm);
    */
	  
    if(list[i]-fatm >= 0)
      {
	if(list[i]-fatm < (latm-fatm+1))
	  {
	    residue->bonded[list[0]-fatm][list[i]-fatm] = nbr[i];
	  }
      }
  }


  return;

}
