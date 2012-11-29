#include "flexaid.h"
#include "boinc.h"

/****************************************************************************/
/****************************************************************************/
/**** This function reads the grid built used and for normal mode       *****/
/**** A vertex in the grid is represented by a combination of all modes *****/
/****************************************************************************/
/****************************************************************************/

void read_normalgrid(FA_Global* FA,char filename[]) {
  FILE* infile_ptr;
  char buffer[150];
  char amplitude[9];

  int i,j;
  
  infile_ptr = NULL;
  if (!OpenFile_B(filename,"r",&infile_ptr))
    Terminate(8);

  FA->normal_grid_points=0;

  
  FA->normal_grid = (float**)malloc(FA->MIN_NORMAL_GRID_POINTS*sizeof(float*));
  if(!FA->normal_grid){
	  fprintf(stderr,"ERROR: memory allocation error for normal_grid\n");
	  Terminate(2);
  }
  memset(FA->normal_grid,NULL,FA->MIN_NORMAL_GRID_POINTS*sizeof(float*));

  for(i=0;i<FA->MIN_NORMAL_GRID_POINTS;i++){
	  FA->normal_grid[i] = (float*)malloc(FA->normal_modes*sizeof(float));
	  if(!FA->normal_grid[i]){
		  fprintf(stderr,"ERROR: memory allocation error for normal_grid[%d].\n", i);
		  Terminate(2);
	  }
	  memset(FA->normal_grid[i],0,FA->normal_modes*sizeof(float));
  }
  
  
  while(fgets(buffer,sizeof(buffer),infile_ptr)!=NULL){
	  // a line represents an intersection of the grid (vertex)
	  // line format: %7.3f %7.3f..'\n'
	  
	  if(FA->normal_grid_points==FA->MIN_NORMAL_GRID_POINTS){
		  FA->MIN_NORMAL_GRID_POINTS*=2;
		  
		  //printf("re-allocating memory for normal_grid\n");
		  FA->normal_grid = (float**)realloc(FA->normal_grid,FA->MIN_NORMAL_GRID_POINTS*sizeof(float*));
		  if(!FA->normal_grid){
			  fprintf(stderr,"ERROR: memory re-allocation error for normal_grid\n");
			  Terminate(2);
		  }
		  memset(&FA->normal_grid[FA->MIN_NORMAL_GRID_POINTS/2],NULL,FA->MIN_NORMAL_GRID_POINTS/2*sizeof(float*));
		  
		  for(i=FA->MIN_NORMAL_GRID_POINTS/2;i<FA->MIN_NORMAL_GRID_POINTS;i++){
			  FA->normal_grid[i] = (float*)malloc(FA->normal_modes*sizeof(float));
			  if(!FA->normal_grid[i]){
				  fprintf(stderr, "ERROR: memory allocation error for normal_grid[%d].\n", i);
				  Terminate(2);
			  }
			  memset(FA->normal_grid[i],0,FA->normal_modes*sizeof(float));
		  }
	  }
	  
	  for(i=0;i<FA->normal_modes;i++){
		  for(j=0;j<8;j++){
			  amplitude[j]=buffer[i*9+j];
		  }
		  amplitude[8]='\0';
		  sscanf(amplitude,"%f",&FA->normal_grid[FA->normal_grid_points][i]);
		  //printf("read col[%d] val[%.4f]\n",i,FA->normal_grid[FA->normal_grid_points][i]);
	  }
	  
	  FA->normal_grid_points++;
  }
  
  //getchar();
  
  CloseFile_B(&infile_ptr,"r");
  
  //printf("normalgrid points=%d\n",FA->normal_grid_points);

  return;
}
