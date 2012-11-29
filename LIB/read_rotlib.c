#include "flexaid.h"
#include "boinc.h"

/***********************************************************************************
 * This subroutine reads a rotamer library and records in an array of structures of
 * the following type:
 *
 * struct RotLib_struct{
 *   char  res[4];
 *   int   nid;
 *   char  name[9];
 *   float pro;
 *   int   nchi;
 *   float chi[4];
 *  };
 ***********************************************************************************/

void read_rotlib(FA_Global* FA,rot** rotamer,char* libfile){
  
  FILE  *infile_ptr;        /* pointer to input file */
  char  buffer[150];        /* a line from the INPUT file */
  char  field[4];           /* field names on INPUT file */
  char  help[20];
  int   i,j;
  //int   k,l,m;
  int   in,dr;              //double range and input starting byte
  char  res[4];
  //int   nid;
  //char  name[9];
  int   obs,total;
  //int   nchi;
  //float chi[4];
  int   ndihedrals = 0;          
  int   rangeflag = 0;          //Phe or Tyr have 2 possible ranges for Chi2

  infile_ptr=NULL;
  if (!OpenFile_B(libfile,"r",&infile_ptr))
    Terminate(8);
  
  dr=0;
  FA->rotlibsize=0;

  while(fgets(buffer, sizeof(buffer),infile_ptr)){
    //printf("%s",buffer);
    for(i=0;i<3;i++){field[i]=buffer[i];}
    field[3]='\0';
    if(strcmp(field,"RES") == 0){
      rangeflag = 0;

      for(i=0;i<3;i++){res[i]=buffer[i+4];}
      res[3]='\0';

      //printf("RESIDUE %s\n",res);
      //ndihedrals=4;
      ndihedrals=number_of_dihedrals(res);
      for(i=0;i<4;i++){help[i]=buffer[i+9];}
      help[4]='\0';
      sscanf(help,"%d",&total);
    }
    if(strcmp(field,"ROT") == 0){
      dr=0;
      (*rotamer)[FA->rotlibsize].numrng++;
      strcpy((*rotamer)[FA->rotlibsize].res,res);
      
      if (strcmp(res,"PHE")==0 || strcmp(res,"TYR")==0) {
	(*rotamer)[FA->rotlibsize].numrng++;
	rangeflag = 1;
      }

      (*rotamer)[FA->rotlibsize].nchi=ndihedrals;
      for(i=0;i<4;i++){help[i]=buffer[i+4];}
      help[4]='\0';
      sscanf(help,"%d",&(*rotamer)[FA->rotlibsize].nid);

      for(i=0;i<9;i++){help[i]=buffer[i+9];}
      help[8]='\0';
      sscanf(help,"%s",(*rotamer)[FA->rotlibsize].name);

      for(i=0;i<4;i++){help[i]=buffer[i+20];}
      help[4]='\0';
      sscanf(help,"%d",&obs);

      (*rotamer)[FA->rotlibsize].tot = 0;
      (*rotamer)[FA->rotlibsize].obs = obs;
      (*rotamer)[FA->rotlibsize].pro = (float)obs/(float)total;


      for(i=0;i<ndihedrals;i++){
	//READ DIHEDRALS
	for(j=0;j<5;j++){help[j]=buffer[25+5*i+j];}
	help[5]='\0';
	sscanf(help,"%f",&(*rotamer)[FA->rotlibsize].chi[i]);
	//printf("CHI[%d]=%6.1f ",i,(*rotamer)[FA->rotlibsize].chi[i]);

	//READ RANGE VALUES
	in=25+5*ndihedrals+1;
	for (j=0;j<10;j++){help[j]=buffer[in+10*i+j];}
	help[10]='\0';
	sscanf(help,"%f %f",&(*rotamer)[FA->rotlibsize].lowr[i][dr],&(*rotamer)[FA->rotlibsize].hghr[i][dr]);
	//printf("LOW[%d]=%6f HIGH[%d]=%6f\n",i,(*rotamer)[FA->rotlibsize].lowr[i][dr],i,(*rotamer)[FA->rotlibsize].hghr[i][dr]);

	if (rangeflag && i==(ndihedrals-1)) { //PHE or TYR (double ranges)
	  dr=1;
	  for (j=0;j<10;j++){help[j]=buffer[in+20*i+j];}
	  help[10]='\0';
	  sscanf(help,"%f %f",&(*rotamer)[FA->rotlibsize].lowr[i][dr],&(*rotamer)[FA->rotlibsize].hghr[i][dr]);
	  //printf("LOW[%d]=%6f HIGH[%d]=%6f\n",i,(*rotamer)[FA->rotlibsize].lowr[i][dr],i,(*rotamer)[FA->rotlibsize].hghr[i][dr]);
	}

	//READ GAUSSIAN BOUNDARIES
	in=in+10*ndihedrals+1+dr*10;
	for (j=0;j<10;j++){help[j]=buffer[in+10*i+j];}
	help[10]='\0';
	//gaus 0 is the positive boundary and 1 the negative boundary
	sscanf(help,"%f %f",&(*rotamer)[FA->rotlibsize].gaus[i][0],&(*rotamer)[FA->rotlibsize].gaus[i][1]);
	if ((*rotamer)[FA->rotlibsize].gaus[i][1]==0){(*rotamer)[FA->rotlibsize].gaus[i][1]=(*rotamer)[FA->rotlibsize].gaus[i][0];}
	//printf("LOW_GAUS[%d]=%6f\tHIGH_GAUS[%d]=%6f\n",i,(*rotamer)[FA->rotlibsize].gaus[i][0],i,(*rotamer)[FA->rotlibsize].gaus[i][1]);
      }

      for(i=0;i<=FA->rotlibsize;i++) {
	if(strcmp((*rotamer)[i].res,res)==0){
	  (*rotamer)[i].tot+=obs;
	  break;
	}
      }
      
      FA->rotlibsize++;
    }
    //PAUSE;
  }  
  CloseFile_B(&infile_ptr,"r");

  return;
}
