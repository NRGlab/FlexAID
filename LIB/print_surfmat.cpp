#include "boinc.h"
#include "flexaid.h"

void print_surfmat(FA_Global* FA,float matrix[9][9], char* filename)
{
  
  FILE* outfile_ptr;
  int i,j;
  const char* ROMAN[8] = {"I","II","III","IV","V","VI","VII","VIII"};
 
  
  if (!OpenFile_B(filename,"w",&outfile_ptr)){
    fprintf(stderr,"ERROR: Could not write to output file: %s\n", filename);
    Terminate(8);
  } 
  
  // header line
  fprintf(outfile_ptr,"          ");
  for(i=1;i<=FA->ntypes;++i) fprintf(outfile_ptr,"%10s",ROMAN[i-1]);
  fprintf(outfile_ptr,"\n");

  for(i=1;i<=FA->ntypes;++i){
    fprintf(outfile_ptr,"%10s",ROMAN[i-1]);
    for(j=1;j<=FA->ntypes;++j){
      fprintf(outfile_ptr,"%10.2f",matrix[i][j]);
    }
    fprintf(outfile_ptr,"\n");
  }
  
  CloseFile_B(&outfile_ptr,"w");

  return;

}
