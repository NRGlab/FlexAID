#include "flexaid.h"
#include "boinc.h"

void read_constraints(FA_Global* FA, atom* atoms, resid* residue, char* filename) {
  
  int i,j;
  char buffer[100];
  char bufnul[12];
  char type[12];
  char *pch;
  char tmp[12];

  int itype=0;
  FILE* infile_ptr;

  
  FA->constraints = (constraint*)malloc(FA->MIN_CONSTRAINTS*sizeof(constraint));
  if(!FA->constraints){
    fprintf(stderr,"ERROR: Could not allocate memory for constraints.\n");
    Terminate(2);
  }

  if (!OpenFile_B(filename,"r",&infile_ptr))
    Terminate(8);
 
  //printf("reading constraint file %s...\n",filename);

  while(fgets(buffer,sizeof(buffer),infile_ptr) != NULL){

    for(i=0;i<8;i++){type[i]=buffer[i];}
    type[8]='\0';

    if(!strcmp(type,"COVALENT")){
      itype=1;
    }else if(!strcmp(type,"INTERACT")){
      itype=0;
    }else{
      continue;
    }


    if(FA->num_constraints == FA->MIN_CONSTRAINTS){
      FA->MIN_CONSTRAINTS += 3;
      
      FA->constraints = (constraint*)realloc(FA->constraints,FA->MIN_CONSTRAINTS*sizeof(constraint));
      if(!FA->constraints){
	fprintf(stderr,"ERROR: Could not re-allocate memory for constraints.\n");
	Terminate(2);
      } 
    }

    FA->constraints[FA->num_constraints].type = itype;
    
    /*
      0         1         2         3         4         5         
      012345678901234567890123456789012345678901234567890123456789
      COVALENT  300 A BTN   904:  88 A SER   558:1.500
    */

    for(i=0;i<4;i++){bufnul[i]=buffer[9+i];}
    bufnul[4]='\0';
    sscanf(bufnul,"%d",&FA->constraints[FA->num_constraints].rnum1);

    for(i=0;i<4;i++){bufnul[i]=buffer[26+i];}
    bufnul[4]='\0';
    sscanf(bufnul,"%d",&FA->constraints[FA->num_constraints].rnum2);
    

    FA->constraints[FA->num_constraints].chn1 = buffer[14] == '-' ? ' ' : buffer[14];
    FA->constraints[FA->num_constraints].chn2 = buffer[31] == '-' ? ' ' : buffer[31];


    for(i=0;i<3;i++){FA->constraints[FA->num_constraints].rnam1[i]=buffer[16+i];}
    FA->constraints[FA->num_constraints].rnam1[3]='\0';

    for(i=0;i<3;i++){FA->constraints[FA->num_constraints].rnam2[i]=buffer[33+i];}
    FA->constraints[FA->num_constraints].rnam2[3]='\0';
    

    for(i=0;i<5;i++){bufnul[i]=buffer[20+i];}
    bufnul[5]='\0';
    sscanf(bufnul,"%d",&FA->constraints[FA->num_constraints].anum1);

    for(i=0;i<5;i++){bufnul[i]=buffer[37+i];}
    bufnul[5]='\0';
    sscanf(bufnul,"%d",&FA->constraints[FA->num_constraints].anum2);


    FA->constraints[FA->num_constraints].bond_len=0.0f;
    if(FA->constraints[FA->num_constraints].type == 1){
      FA->constraints[FA->num_constraints].bond_len = (float)atof(&buffer[43]);

    }else{
	    
            sscanf(&buffer[43],"%s",tmp);
	    
	    FA->constraints[FA->num_constraints].force_interaction = 0;
	    if((pch=strchr(tmp,'f')) != NULL){
		    FA->constraints[FA->num_constraints].force_interaction = 1;
		    *pch=' ';
	    }
	    sscanf(tmp,"%f",&FA->constraints[FA->num_constraints].interaction_factor);
	    
	    /*
	    printf("New interaction factor: %.2f with force %s\n", 
		   FA->constraints[FA->num_constraints].interaction_factor, 
		   FA->constraints[FA->num_constraints].force_interaction?"ON":"OFF");
	    */
    }
	    
    FA->constraints[FA->num_constraints].id = FA->num_constraints;

	    
    if(assign_constraint(FA,atoms,residue,&FA->constraints[FA->num_constraints])){
	    ++FA->num_constraints;
    }

  }

  
  for(j=0;j<FA->num_constraints;++j){
	  update_constraint(atoms,FA->constraints[j].inum1,&FA->constraints[j]);
	  update_constraint(atoms,FA->constraints[j].inum2,&FA->constraints[j]);
  }


  CloseFile_B(&infile_ptr,"r");

}
