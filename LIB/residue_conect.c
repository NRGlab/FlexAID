#include "flexaid.h"
#include "boinc.h"

/******************************************************************************
 * SUBROUTINE residue_conect reads and creates conectivity for protein residues
 * it also reads the flexible dihedral angles definitions.
 ******************************************************************************/
void residue_conect(FA_Global* FA,atom* atoms,resid* residue,char aminofile[]){

	FILE *infile_ptr;        /* pointer to input file */
	char buffer[81];         /* a line from the INPUT file */
	char bufnul[7];          /* dumb string to read input strings */
	char rnam[4];            /* residue name */
	//char rnamrev[4];

	int  i,j,l;              /* dumb counter */
	int  flag;               /* flag for reading conect data */
	//int  type;               /* temporary atom type */
	int  natm=0;               /* counter of number of atoms in ideal residue */
	int  natm_tmp;
	char field[7];
	char num_char[6];        /* string to read integer */
	int atm;                 /* number of the atom we read info */
	int atn;                 /* atn is covelently bonded to atm */
	int conect[25][5];       /* conectivity matrix of ideal residue */
	int dihed[10];           /* conectivity matrix of ideal residue */
	int fdih=0;
	int dihflag;

	infile_ptr=NULL;
	if (!OpenFile_B(aminofile,"r",&infile_ptr))
		Terminate(8);

	flag=0;

	for(j=0;j<=9;j++){
		dihed[j]=0;
	}
	while (fgets(buffer, sizeof(buffer),infile_ptr)){
		for(i=0;i<=5;i++){field[i]=buffer[i];}
		field[6]='\0';

		if(strcmp(field,"RESIDU") == 0 || strcmp(field,"NUCLEO") == 0){ 

			flag=0;
			natm=0;
			fdih=0;

			for(i=0;i<=2;i++){bufnul[i]=buffer[i+7];}
			bufnul[3]='\0';
			strcpy(rnam,bufnul);
			//printf("%s\n",rnam);
			
			/* RNA modifs.
			rnamrev[0]=rnam[2];
			rnamrev[1]=rnam[1];
			rnamrev[2]=rnam[0];
			rnamrev[3]='\0';
			*/
		}
    

		if(strcmp(field,"CONECT") == 0){
			natm++;
			conect[natm][0]=0;
			for(i=0;i<=1;i++){num_char[i]=buffer[7+i];}
			num_char[2]='\0';
			sscanf(num_char,"%d",&atm);
			for(j=0;j<=3;j++){
				atn=0;
				for(i=0;i<=1;i++){num_char[i]=buffer[10+j*3+i];}
				num_char[2]='\0';
				sscanf(num_char,"%d",&atn);
				if(atn != 0){
					conect[natm][j+1]=atn;
					conect[natm][0]++;
				}else{
					break;
				}
			}
		}

		if(strcmp(field,"FLEDIH") == 0){
			fdih++;
			for(i=0;i<=1;i++){
				num_char[i]=buffer[10+i];
			}
			num_char[2]='\0';
			sscanf(num_char,"%d",&dihed[fdih]);
			/*printf("fledih:%d %d or %s\n",fdih,dihed[fdih],num_char);
			  PAUSE;*/
			dihflag=1;
		}


		if(strcmp(field,"RESEND") == 0){ 
			flag=1;
		}

  
		if(flag == 1){    
			for(i=1;i<=FA->res_cnt;i++){
				natm_tmp = natm;
				
				if(strcmp(residue[i].name,rnam)==0){ // RNAmodifs: || (!FA->is_protein && strcmp(residue[i].name,rnamrev)==0)){
					//printf("natm_tmp[res=%d]=[%d]\n",residue[i].number,natm_tmp);

					for(j=0;j<FA->MIN_FLEX_BONDS;j++)
						residue[i].bond[j]=0;
	  
					if(residue[i].latm[0]-residue[i].fatm[0]+1 == natm_tmp){
						/*
						printf("Residue name: %s number: %d ok!\n", residue[i].name, i);
						printf("First atom: %d(%s) Last atom: %d(%s)\n", residue[i].fatm[0], atoms[residue[i].fatm[0]].name, residue[i].latm[0], atoms[residue[i].latm[0]].name);
						*/

						for(l=1;l<=natm_tmp;l++){
							atoms[residue[i].fatm[0]+l-1].bond[0]=conect[l][0];
							for(j=1;j<=atoms[residue[i].fatm[0]+l-1].bond[0];j++){
								atoms[residue[i].fatm[0]+l-1].bond[j]=conect[l][j]+residue[i].fatm[0]-1;
							}
						}

						if(i != 1){
							if(residue[i-1].type == 0){
								j=residue[i-1].fatm[0];
								while(j<=residue[i-1].latm[0]){
									if(strcmp(atoms[j].name," C  ") == 0){
										atoms[residue[i].fatm[0]].bond[atoms[residue[i].fatm[0]].bond[0]+1]=j;
										atoms[residue[i].fatm[0]].bond[0]++;
										break;
									}
									j++;
								}
							}
						}
						if(i != FA->res_cnt) {
							if(residue[i+1].type == 0){
								j=residue[i+1].fatm[0];
								while(j<=residue[i+1].latm[0]){
									if(strcmp(atoms[j].name," N  ") == 0){
										atoms[residue[i].fatm[0]+2].bond[atoms[residue[i].fatm[0]+2].bond[0]+1]=j;
										atoms[residue[i].fatm[0]+2].bond[0]++;
										break;
									}
									j++;
								}
							}
						}

						residue[i].fdih=fdih;
						for(j=1;j<=fdih;j++){
							residue[i].bond[j] = residue[i].fatm[0] + dihed[j] - 1;
	      
							/*
							  if(i==2){
							  printf("fatm[0]=[%d]\tdihed[%d]=%d\tresidue[%d].bond[%d]=%d\n",residue[i].fatm[0],j,dihed[j],i,j,residue[i].bond[j]);
							  PAUSE;
							  }
							*/
						}
					}else{
						printf("The residue %s%d%c misses atoms: FATM=%d\tLATM=%d\n",residue[i].name,
						       residue[i].number,residue[i].chn,residue[i].fatm[0],residue[i].latm[0]);
						//PAUSE;
						//exit(8);
					}
				}
			}
			flag=0;
		}
	}
	CloseFile_B(&infile_ptr,"r");

	return;
}
