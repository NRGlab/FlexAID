#include "flexaid.h"
#include "boinc.h"

////////////////////////////////////////////////
//// reads fixed conformation derived from  ////
//// the Holo form of the HAP2 database     ////
////////////////////////////////////////////////

void read_rotobs(FA_Global* FA,rot** rotamer,char* filename)
{
  int i,j;
  FILE* infile_ptr;
  char buffer[MAX_PATH__*2];
  char residue[4];
  char chi[9];
  int nchi;
  int nid;

  infile_ptr=NULL;
  if (!OpenFile_B(filename,"r",&infile_ptr))
    Terminate(8);

  FA->rotlibsize=0;
  nid=1;

  while(fgets(buffer,sizeof(buffer),infile_ptr)!=NULL){
    for(i=0;i<3;i++){residue[i]=buffer[i];}
    residue[3]='\0';

    nchi=number_of_dihedrals(residue);
    
    // residue is not a flexible residue ( no CHI angle(s) )
    if(nchi==0)
      continue;

    strcpy((*rotamer)[FA->rotlibsize].res,residue);
    for(i=0;i<nchi;i++){
      for(j=0;j<8;j++){
	chi[j]=buffer[j+i*9+4];
      }
      chi[8]='\0';

      sscanf(chi,"%f",&(*rotamer)[FA->rotlibsize].chi[i]);
      (*rotamer)[FA->rotlibsize].obs=1;
      (*rotamer)[FA->rotlibsize].tot=0;
    }
    
    for(i=0;i<=FA->rotlibsize;i++){
      if(strcmp((*rotamer)[i].res,residue)==0){
	(*rotamer)[i].tot++;
	(*rotamer)[FA->rotlibsize].nid=(*rotamer)[i].tot;
	break;
      }
    }
    

    FA->rotlibsize++;
    
    if(FA->rotlibsize == FA->MIN_ROTAMER_LIBRARY_SIZE){
      //printf("reallocating memory for rotamer\n");

      FA->MIN_ROTAMER_LIBRARY_SIZE += 100;

      (*rotamer) = (rot*)realloc((*rotamer),FA->MIN_ROTAMER_LIBRARY_SIZE*sizeof(rot));
      if(!(*rotamer)){
	fprintf(stderr,"ERROR: memory allocation error for rotamer.\n");
	Terminate(2);
      }

      memset(&(*rotamer)[FA->MIN_ROTAMER_LIBRARY_SIZE-100],0,100*sizeof(rotamer));

      //printf("memory reallocated for rotamer\n");
    }
  }

  CloseFile_B(&infile_ptr,"r");

  return;
}
