#include "flexaid.h"
#include "boinc.h"
/***************************************************************************** 
 * SUBROUTINE read coor, gets a coordinates line from read_pdb and extracts
 * the coordinates, atom name, chain name, counts the number of ligands, 
 * residues, the initial and final atom of each residue and hetero group
 *****************************************************************************/
void read_coor(FA_Global* FA,atom** atoms,resid** residue,char line[], char res_numold[]){
	char  name[7];             /* 6 letter code of field name on PDB file, e.g. HETATM */
	char  coor_char[10];        /* string used to read the coordinates */
	char  num_char[6];         /* string to read the atom number field */
	char  atm_typ[5];          // temporary type
	char  res_new[4];          // temporary residue name 
	char  res_num[5];          // temporary residue number from the PDB annotation 
	//char  res_numold[5];       // another temporary residue number from the annot.
  
	int   i,j;                  /* dumb counters */

	/*

	  01234567890123456789012345678901234567890123456789012345678901234567890123456789
	  ATOM    239  CB  ALA A  46      17.761  -3.260 -12.974  1.00 20.41           C  
    
	*/
  
	for(i=0;i<=5;i++){name[i]=line[i];}
	name[6]='\0';
  
	for(i=0;i<=3;i++){atm_typ[i]=line[i+12];}
	atm_typ[4]='\0';

	if(line[16]==' ' || line[16]=='A'){
		FA->atm_cnt++;
    
		if(FA->atm_cnt==FA->MIN_NUM_ATOM){
			//printf("re-allocating memory for atoms\n");
			FA->MIN_NUM_ATOM*=2;

			(*atoms) = (atom*)realloc((*atoms),FA->MIN_NUM_ATOM*sizeof(atom));
      
			if(!(*atoms)){
				fprintf(stderr,"ERROR: memory allocation error for atoms.\n");
				Terminate(2);
			}
      
			memset(&(*atoms)[FA->MIN_NUM_ATOM/2],0,FA->MIN_NUM_ATOM/2*sizeof(atom));
			//printf("memory re-allocated for atoms\n");
		}
    
		(*atoms)[FA->atm_cnt].eigen = NULL;
		(*atoms)[FA->atm_cnt].ncons=0;
		(*atoms)[FA->atm_cnt].cons=NULL;
		(*atoms)[FA->atm_cnt].optres=NULL;
		(*atoms)[FA->atm_cnt].graph=0;
		(*atoms)[FA->atm_cnt].coor_ref=NULL;
        
		strcpy((*atoms)[FA->atm_cnt].name,atm_typ);
		if(strcmp((*atoms)[FA->atm_cnt].name," OXT")==0){
			(*residue)[FA->res_cnt].ter = 1;
			//printf("Residue Ter: %d\n", (*residue)[FA->res_cnt].ter);
		}

		(*atoms)[FA->atm_cnt].isbb=0;

		if (!strcmp((*atoms)[FA->atm_cnt].name," CB ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," CA ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," N  ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," O  ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," C  ") ||
		    !strcmp((*atoms)[FA->atm_cnt].name," OXT")) 
		{ (*atoms)[FA->atm_cnt].isbb=1; }
    

		(*atoms)[FA->atm_cnt].radius=assign_radius((*atoms)[FA->atm_cnt].name);
    
        (*atoms)[FA->atm_cnt].type = atoi(&line[76]);
        
		for(j=0;j<=4;j++){num_char[j]=line[j+6];}
		num_char[5]='\0';
		sscanf(num_char,"%d",&i);
    
		(*atoms)[FA->atm_cnt].number=i;

		/* maps the PDB num into the internal counter */
		FA->num_atm[i]=FA->atm_cnt;

		for (j=0;j<=2;j++){
			for(i=0;i<=7;i++){
				coor_char[i]=line[30+i+j*8];
			}
			coor_char[8]='\0';
			sscanf(coor_char,"%f",&(*atoms)[FA->atm_cnt].coor[j]);
			sscanf(coor_char,"%f",&(*atoms)[FA->atm_cnt].coor_ori[j]);
		}
    
		for(i=0;i<=2;i++){res_new[i]=line[i+17];}
		res_new[3]='\0';
		for(i=0;i<=3;i++){res_num[i]=line[i+22];}
		res_num[4]='\0';
    
		if(strcmp(res_new,(*residue)[FA->res_cnt].name) != 0   /* change of res name        */
		   || strcmp(res_num,res_numold) != 0           /* change of res number      */
			){
			FA->res_cnt++;

			if(FA->res_cnt==FA->MIN_NUM_RESIDUE){
				//printf("re-allocating memory for residue\n");
				FA->MIN_NUM_RESIDUE *= 2;
				(*residue) = (resid*)realloc((*residue),FA->MIN_NUM_RESIDUE*sizeof(resid));
				if(!(*residue)){
					fprintf(stderr,"ERROR: memory allocation error for residue.\n");
					Terminate(2);
				}
				memset(&(*residue)[FA->MIN_NUM_RESIDUE/2],0,FA->MIN_NUM_RESIDUE/2*sizeof(residue));
				//printf("memory re-allocated for residue\n");
			}

			(*residue)[FA->res_cnt].fatm = (int*)malloc(FA->MIN_ROTAMER*sizeof(int));
			(*residue)[FA->res_cnt].latm = (int*)malloc(FA->MIN_ROTAMER*sizeof(int));
			(*residue)[FA->res_cnt].bond = (int*)malloc(FA->MIN_FLEX_BONDS*sizeof(int));
			if(!(*residue)[FA->res_cnt].fatm ||
			   !(*residue)[FA->res_cnt].latm ||
			   !(*residue)[FA->res_cnt].bond){
				fprintf(stderr,"ERROR: memory allocation error for residue.fatm || .latm || .bond.\n");
				Terminate(2);
			}

			memset((*residue)[FA->res_cnt].fatm,0,FA->MIN_ROTAMER*sizeof(int));
			memset((*residue)[FA->res_cnt].latm,0,FA->MIN_ROTAMER*sizeof(int));
			memset((*residue)[FA->res_cnt].bond,0,FA->MIN_ROTAMER*sizeof(int));

			(*residue)[FA->res_cnt].ter = 0;
			(*residue)[FA->res_cnt].rot=0;
			(*residue)[FA->res_cnt].fdih=0;
			(*residue)[FA->res_cnt].trot = 0;
			(*residue)[FA->res_cnt].gpa=NULL;
			(*residue)[FA->res_cnt].bonded=NULL;
			(*residue)[FA->res_cnt].fatm[0]=FA->atm_cnt;     
			(*residue)[FA->res_cnt-1].latm[0]=FA->atm_cnt-1;
			(*residue)[FA->res_cnt].chn=line[21];
			(*residue)[FA->res_cnt].ins=line[26];
      
			strcpy((*residue)[FA->res_cnt].name,res_new);
			strcpy(res_numold,res_num);

			if(strcmp(name,"ATOM  ")==0) {
				(*residue)[FA->res_cnt].type=0;                       /* protein residue */
			}else if(strcmp(name,"HETATM")==0) {
				(*residue)[FA->res_cnt].type=1;
			}

			sscanf(res_num,"%d",&(*residue)[FA->res_cnt].number);

			//printf("New residue: %d fatm[%d]=%d(%s-%s) latm[%d]=%d(%s-%s) :: %d\n", FA->res_cnt, FA->res_cnt, (*atoms)[(*residue)[FA->res_cnt].fatm[0]].number, (*residue)[FA->res_cnt].name, (*atoms)[(*residue)[FA->res_cnt].fatm[0]].name, FA->res_cnt-1, (*atoms)[(*residue)[FA->res_cnt-1].latm[0]].number, (*residue)[FA->res_cnt-1].name, (*atoms)[(*residue)[FA->res_cnt-1].latm[0]].name,(*residue)[FA->res_cnt].number);
            
		}

		(*atoms)[FA->atm_cnt].ofres=FA->res_cnt;

	}

	return;
}
