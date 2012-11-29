#include "flexaid.h"
#include "boinc.h"

/******************************************************************************
 * SUBROUTINE assign_types assign atom types for protein and hetero group atoms
 ******************************************************************************/
void assign_types(FA_Global* FA,atom* atoms,resid* residue,char aminofile[]){

	FILE *infile_ptr;        /* pointer to input file */
	char buffer[81];         /* a line from the INPUT file */
	char bufnul[7];          /* dumb string to read input strings */
	int  i,j,k;                /* dumb counter */
	char anam[5];            /* temporary atom name */
	char rnam[4],rnamrev[4]; /* temporary residue name */
	int  type;               /* temporary atom type */
	//char field[7];
	char recs;
	int  rec[4];
	int rot;

	infile_ptr=NULL;
	if (!OpenFile_B(aminofile,"r",&infile_ptr))
		Terminate(8);

	while (fgets(buffer, sizeof(buffer),infile_ptr)){

		for(i=0;i<=5;i++){bufnul[i]=buffer[i];}
		bufnul[6]='\0';

		if(strcmp(bufnul,"RESIDU") == 0 || strcmp(bufnul,"NUCLEO") == 0){ 
			for(i=0;i<=2;i++){rnam[i]=buffer[i+7];}
			rnam[3]='\0';

			rnamrev[0]=rnam[2];
			rnamrev[1]=rnam[1];
			rnamrev[2]=rnam[0];
			rnamrev[3]='\0';
			//printf("%s\n",rnam);
		}

		if(strcmp(bufnul,"ATMTYP") == 0){ 
			bufnul[0]=buffer[10];
			bufnul[1]=buffer[11];
			bufnul[2]='\0';
			sscanf(bufnul,"%d",&type);
      
			for(i=0;i<=3;i++){anam[i]=buffer[i+12];}
			anam[4]='\0';
      
			if(type > FA->ntypes){
				printf("WARNING: res %s atom %s has atom type %d when %d types are defined\n",
				       rnam, anam, type, FA->ntypes);
				printf("WARNING: type %d is set to neutral (6)\n", type);

				type = 6;
			}


			recs=buffer[17];

			/*printf("recs:%c\n",recs);*/

			if(recs == 'm'){
				for(i=0;i<=3;i++){
					for(j=0;j<=2;j++){bufnul[j]=buffer[18+i*3+j];}
					bufnul[3]='\0';
					sscanf(bufnul,"%d",&rec[i]);
				}
				/*printf("%d %d %d\n",rec[0],rec[1],rec[2]);*/
			}
      
			/*PAUSE;*/
			for(k=1;k<=FA->res_cnt;k++){
				rot=residue[k].rot;
				for(i=residue[k].fatm[rot];i<=residue[k].latm[rot];i++){

					if(
						(strcmp(residue[atoms[i].ofres].name,rnam) == 0 ||
						 (!FA->is_protein && strcmp(residue[atoms[i].ofres].name,rnamrev) == 0)) &&
						strcmp(atoms[i].name,anam) == 0){
						
						atoms[i].type=type;
						atoms[i].recs=recs;
						
						if(atoms[i].recs == 'm'){
							for(j=0;j<=2;j++){
								atoms[i].rec[j]=rec[j]+
									residue[atoms[i].ofres].fatm[residue[atoms[i].ofres].rot]-1;
							}
							if(rec[3] != 0){
								atoms[i].rec[3]=rec[3]+
									residue[atoms[i].ofres].fatm[residue[atoms[i].ofres].rot]-1;
							}else{
								atoms[i].rec[3]=0;
							}
						}else{
							for(j=0;j<=3;j++){
								atoms[i].rec[j]=0;
							}
						}
					}
				}
			}
		}
	}
	CloseFile_B(&infile_ptr,"r");

	// Set a Hydrophilic type to water molecules
	for(k=1;k<=FA->res_cnt;k++){ 
		rot=residue[k].rot;
		for(i=residue[k].fatm[rot];i<=residue[k].latm[rot];i++){ 
			if(strcmp(residue[atoms[i].ofres].name,"HOH") == 0){atoms[i].type = 1;}
		}
	}

	// Set a '' type for Metals in general
	for(k=1;k<=FA->res_cnt;k++){ 
		rot=residue[k].rot;
	   
		for(i=residue[k].fatm[rot];i<=residue[k].latm[rot];i++){ 
			if(strcmp(residue[atoms[i].ofres].name,"FE ") == 0 ||
			   strcmp(residue[atoms[i].ofres].name,"CU ") == 0 ||
			   strcmp(residue[atoms[i].ofres].name,"ZN ") == 0 ||
			   strcmp(residue[atoms[i].ofres].name,"MG ") == 0 ||
			   strcmp(residue[atoms[i].ofres].name,"CO ") == 0 ||
			   strcmp(residue[atoms[i].ofres].name,"NI ") == 0 ||
			   strcmp(residue[atoms[i].ofres].name,"MN ") == 0){
				atoms[i].type = FA->metaltype;

				printf("metal %s%d%c set to type %d\n",
				       residue[atoms[i].ofres].name,
				       residue[atoms[i].ofres].number,
				       residue[atoms[i].ofres].chn,
				       FA->metaltype);

			}
		}
	}

	return;
}
