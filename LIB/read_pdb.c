#include "flexaid.h"
#include "boinc.h"

/******************************************************************************
 * SUBROUTINE read_pdb opens a pdb file, reads is and sends each ATOM or  
 * HETATM line to the subroutine read_coor that will read the information 
 * on each line                                                           
 *****************************************************************************/
                         
void read_pdb(FA_Global* FA,atom** atoms,resid** residue, char* pdb_name){
  
  FILE *infile_ptr;        /* pointer to PDB file */
  char buffer[85];         /* input line from PDB file */
  char name[7];            /* 6 letter field name, e.g. HETATM or ATOM */
  char resname[4];         /* residue name */
  int i;                   /* dumb counter */
  char resnumold[5];

  // dynamic memory allocation for atoms and residue
  (*atoms) = (atom*)malloc(FA->MIN_NUM_ATOM*sizeof(atom));
  (*residue) = (resid*)malloc(FA->MIN_NUM_RESIDUE*sizeof(resid));
  FA->num_atm = (int*)malloc(100000*sizeof(int));

  if(!(*atoms) || !(*residue) || !FA->num_atm){
    fprintf(stderr, "ERROR: memory allocation error for atoms || residue.\n");
    Terminate(2);
  }
  
  memset((*atoms),0,FA->MIN_NUM_ATOM*sizeof(atom));
  memset((*residue),0,FA->MIN_NUM_ATOM*sizeof(residue));
  memset(FA->num_atm,0,100000*sizeof(int));
  
  (*atoms)[0].eigen = NULL;
  (*atoms)[0].cons = NULL;
  (*atoms)[0].optres = NULL;
  (*atoms)[0].par = NULL;
  (*residue)[0].gpa = NULL;
  
  (*residue)[0].fatm = (int*)malloc(sizeof(int));
  (*residue)[0].latm = (int*)malloc(sizeof(int));
  if(!(*residue)[0].fatm || !(*residue)[0].latm){
	  fprintf(stderr,"ERROR: Could not allocate memory for residue[0].fatm\n");
	  Terminate(2);
  }
  
  (*residue)[0].bond = NULL;
  (*residue)[0].bonded = NULL;

  strcpy((*residue)[0].name,"   ");
  
  infile_ptr=NULL;
  
  infile_ptr = fopen(pdb_name,"r");
  if(infile_ptr == NULL){
	  fprintf(stderr,"ERROR: Could not read temporary PDB file %s\n", pdb_name);
	  Terminate(8);
  }

  /*
  if (!OpenFile_B(pdb_name,"r",&infile_ptr))
    Terminate(8);
  */

  strcpy(resnumold,"-999");
  while (fgets(buffer, sizeof(buffer),infile_ptr)!=NULL){
    for(i=0;i<=5;i++){name[i]=buffer[i];}
    name[6]='\0';

    for(i=0;i<=2;i++){resname[i]=buffer[17+i];}
    resname[3]='\0';

    if( strcmp(name,"ATOM  ") == 0){
      read_coor(FA,atoms,residue,buffer,resnumold);
    }else if(strcmp(name,"HETATM") == 0){
      read_coor(FA,atoms,residue,buffer,resnumold);
    }else if(strcmp(name,"CONECT") == 0){
      read_conect(FA,atoms,buffer);
    }
  }

  fclose(infile_ptr);
  
  //CloseFile_B(&infile_ptr,"r");
  
  return;
}
