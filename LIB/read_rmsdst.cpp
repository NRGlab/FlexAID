#include "flexaid.h"
#include "boinc.h"

/*****************************************************************************
 * SUBROUTINE read_lig reads the ligand information from the input file 
 *****************************************************************************/
int read_rmsdst(FA_Global* FA,atom* atoms,resid* residue,char rmsdst_file[]){

	FILE *infile_ptr;         /* pointer to PDB file */
	char buffer[100];        /* input line from PDB file */
	char field[7];           /* field names on INPUT file */
	
	char rres_[4], rnum_[5], anum_[6], rchn_, coor_[9];
	int rnum, anum;
	int i,j,k,l;
	int rmsd_atoms = 0;

	if (!OpenFile_B(rmsdst_file,"r",&infile_ptr))
		Terminate(8);

	while(fgets(buffer, sizeof(buffer), infile_ptr) != NULL){
		
		strncpy(field, buffer, 6);
		field[6]='\0';

		//0         1         2         3         4         5         6         
		//0123456789012345678901234567890123456789012345678901234567890123456789
		//ATOM    187  CG  LEU A  39      18.853   9.261   3.766  1.00 20.15           C  
		if(!strcmp(field,"HETATM") || !strcmp(field,"ATOM  ")){
			
			strncpy(anum_,&buffer[6],5); 
			anum_[5]='\0';
			anum = atoi(anum_);

			strncpy(rres_,&buffer[17],3); 
			rres_[3]='\0';
			strncpy(rnum_,&buffer[22],4); 
			rnum_[4]='\0';

			rnum = atoi(rnum_);
			rchn_ = buffer[21];
			
			for(i=1; i<=FA->res_cnt; i++){

				if(!strcmp(residue[i].name,rres_) &&
				   residue[i].number == rnum &&
				   residue[i].chn == rchn_){

					for(j=residue[i].fatm[0]; j<=residue[i].latm[0]; j++){
						if(atoms[j].number == anum){

							atoms[j].coor_ref = (float*)malloc(3*sizeof(float));
							if(atoms[j].coor_ref == NULL){
								fprintf(stderr,"ERROR: Could not allocate memory for coor_ref\n");
								Terminate(2);
							}
							
							for(k=0; k<3; k++){
								for(l=0; l<8; l++){
									coor_[l] = buffer[30+l+k*8];
								}
								coor_[8] = '\0';
								sscanf(coor_, "%f", &atoms[j].coor_ref[k]);
							}
							
							rmsd_atoms++;
						}
					}
				}
			}
		}

	}

	CloseFile_B(&infile_ptr,"r");

	return rmsd_atoms;
}
