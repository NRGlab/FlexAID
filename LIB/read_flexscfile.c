#include "flexaid.h"
#include "boinc.h"

/***********************************************************************************
 * This subroutine reads the identities of the side chains to be set flexible
 * the input file looks like this:
 ***********************************************************************************/

void read_flexscfile(FA_Global* FA,resid* residue,rot** rotamer,char file[], char* rotlib, char* rotobs){

  FILE  *infile_ptr;        /* pointer to input file */
  char buffer[MAX_PATH__*2];         /* a line from the INPUT file */
  char field[7];
  char number[5];
  char residname[4];
  int i,j;    

  // memory allocation for rotamer
  (*rotamer) = (rot*)malloc(FA->MIN_ROTAMER_LIBRARY_SIZE*sizeof(rot));
  if(!(*rotamer)){
    fprintf(stderr,"ERROR: memory allocation error for rotamer\n");
   Terminate(2);
  }
  
  memset((*rotamer),0,FA->MIN_ROTAMER_LIBRARY_SIZE*sizeof(rot));

  // memory allocation for flex_res
  FA->flex_res = (flxsc*)malloc(FA->MIN_FLEX_RESIDUE*sizeof(flxsc));
  if(!FA->flex_res){
    fprintf(stderr,"ERROR: memory allocation error for flex_res\n");
    Terminate(2);
  }
  
  memset(FA->flex_res,0,FA->MIN_FLEX_RESIDUE*sizeof(flxsc));

  ////////////////////////////////////

  infile_ptr=NULL;
  if(!OpenFile_B(file,"r",&infile_ptr))
    Terminate(8);

  FA->nflxsc=0;

  while(fgets(buffer, sizeof(buffer),infile_ptr)){

    for (i=0;i<6;i++){field[i]=buffer[i];}
    field[6]='\0';

    if(strcmp(field,"RESIDU") == 0){
      for(i=14;i<17;i++){
          residname[i-14]=buffer[i];
      }
      residname[3]='\0';
      
      //01234567890123456789012345678901234567890123456789
      //RESIDU  110 A LEU
      if (strcmp(residname,"GLY")==0 || strcmp(residname,"ALA")==0 || strcmp(residname,"PRO")==0) 
	continue;

      if(FA->nflxsc==FA->MIN_FLEX_RESIDUE){
	FA->MIN_FLEX_RESIDUE+=5;
	
	FA->flex_res = (flxsc*)realloc(FA->flex_res,FA->MIN_FLEX_RESIDUE*sizeof(flxsc));
	if(!FA->flex_res){
	  fprintf(stderr,"ERROR: memory allocation error for flex_res\n");
	  Terminate(2);
	}

	memset(&FA->flex_res[FA->MIN_FLEX_RESIDUE-5],0,5*sizeof(flxsc));
      }
  
      strcpy(FA->flex_res[FA->nflxsc].name,residname);
      
      FA->flex_res[FA->nflxsc].chn=buffer[12];
 
      for(i=7;i<11;i++){number[i-7]=buffer[i];}
      number[4]='\0';
      sscanf(number,"%d",&FA->flex_res[FA->nflxsc].num);
      
      FA->flex_res[FA->nflxsc].chn = buffer[12] == '-' ? ' ' : buffer[12];

      FA->nflxsc++;

      //printf("FLEXSC(%d): %4d %c %3s\n",FA->nflxsc,FA->flex_res[FA->nflxsc].num,
      //     FA->flex_res[FA->nflxsc].chn,FA->flex_res[FA->nflxsc].name);


      // also set inum
      
      j=1;
      while(1){
	
	if((strcmp(residue[j].name,FA->flex_res[FA->nflxsc-1].name)==0 &&
	    residue[j].chn == FA->flex_res[FA->nflxsc-1].chn &&
	    residue[j].number == FA->flex_res[FA->nflxsc-1].num) || 
	   j > FA->res_cnt)
	  break;
	
	j++;

      }

      if(j <= FA->res_cnt){
	
	FA->flex_res[FA->nflxsc-1].inum = j;
	
	set_intprob(&FA->flex_res[FA->nflxsc-1]);

	//	printf("Residue: %s %c %4d Found\n",residue[j].name,
	//     residue[j].chn,residue[j].number);
	

      }else{

	FA->nflxsc--;
	
	//printf("Residue: %s %c %4d Not Found\n",FA->flex_res[FA->nflxsc].name,
	//     FA->flex_res[FA->nflxsc].chn,FA->flex_res[FA->nflxsc].num);

      }


    }else if(strcmp(field,"PROBLT") == 0){

      /*
	for (i=7;i<10;i++) {
	residname[i-7]=buffer[i];
	}
	residname[3]='\0';
	
	if (strcmp(residname,"GLY")==0 || strcmp(residname,"ALA")==0 || strcmp(residname,"PRO")) 
	continue;
	
	for (i=11;i<16;i++) {
	probability[i-11]=buffer[i];
	}
	probability[5]='\0';
	
	for (i=0;i<17;i++){
	if(strcmp(FA->flex_prob[i].name,residname)==0){
	sscanf(probability, "%f", &FA->flex_prob[i].prob);
	break;
	}
	}
      */

    }else if(strcmp(field,"ROTOBS")==0){

      //sscanf(buffer, "%s %s",field,rotobs);
      
      // *** NEW in >= 1.27
      FA->rotobs = 1;

    }

  }

  CloseFile_B(&infile_ptr,"r");


  return;
}
