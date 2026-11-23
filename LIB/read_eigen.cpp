#include "flexaid.h"
#include "boinc.h"

/***********************************************************/
/****** This function reads the eigen vectors used   *******/
/****** by the normal mode simulation                *******/
/***********************************************************/

void read_eigen(FA_Global *FA,char filename[]){
  FILE* infile_ptr;
  char buffer[150];
  char vector[9];

  int eigensize;
  int i,j,k;
    
  if (!OpenFile_B(filename,"r",&infile_ptr))
    Terminate(8);

  FA->num_eigen=0;
  FA->normal_nodes=0;

  //////////////////////////////////////
  ///// memory allocation for eigenvector

  eigensize = 3*FA->MIN_NUM_ATOM;

  FA->eigenvector = (float**)malloc(eigensize*sizeof(float*));
  if(!FA->eigenvector){
    fprintf(stderr,"ERROR: memory allocation error for eigenvector\n");
    Terminate(2);
  }
  memset(FA->eigenvector,NULL,eigensize*sizeof(float*));

  for(i=0;i<eigensize;i++){
    FA->eigenvector[i] = (float*)malloc(FA->normal_modes*sizeof(float));
    if(!FA->eigenvector[i]){
      fprintf(stderr,"ERROR: memory allocation error for eigenvector[%d].\n",i);
      Terminate(2);
    }

    memset(FA->eigenvector[i],0,FA->normal_modes*sizeof(float));
  }

  //////////////////////////////////////////

  k=0;
  while(fgets(buffer,sizeof(buffer),infile_ptr)!=NULL){
    // a line represents the vectors for each mode
    // it only considers C,CA and N backbone atoms (super mode)
    // O uses the same eigen than N
    // side-chain uses the same eigen than CA
    
    for(i=0;i<FA->normal_modes;i++){
      for(j=0;j<8;j++){
	vector[j]=buffer[i*9+j];
      }
      vector[8]='\0';
      sscanf(vector,"%f",&FA->eigenvector[FA->num_eigen][i]);
      //printf("read %.5f (num_eigen=%d) normal_modes=%d\n",FA->eigenvector[FA->num_eigen][i],FA->num_eigen,FA->normal_modes);
    }
    
    FA->num_eigen++;

    k++;
    if(k % 3 == 0)
      FA->normal_nodes++;
  }
  
  CloseFile_B(&infile_ptr,"r");

  printf("normal_nodes=%d\n",FA->normal_nodes);
  //getchar();

  return;
}
