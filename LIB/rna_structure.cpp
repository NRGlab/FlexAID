/*
 *  rna_structure.c
 *  FlexAID_Xcode
 *
 *  Created by Francis Gaudreault on 13-01-15.
 *  Copyright 2013 Universite de Sherbrooke. All rights reserved.
 *
 */

#include "flexaid.h"
#include "boinc.h"

/****************************************************/
/******  This function reads the input target *******/
/******  file and if ATOM records belong to   *******/
/******  nucleic acid residues in average     *******/
/******  it is considered as a RNA structure  *******/
/****************************************************/

int rna_structure(char* infile)
{
    FILE* infile_ptr = NULL;
    
    char res[4];
    char buffer[80];
    
    int n_atoms = 0;  // atom records counter
    int n_naa = 0; // nucleic acid atom counter
    int n_aaa = 0; // amino acid atom counter
    
    if(!OpenFile_B(infile,"r",&infile_ptr)){
		fprintf(stderr,"ERROR: Could not read PDB file %s in rna_structure\n", infile);
		Terminate(20);
	}
    
	while(fgets(buffer,sizeof(buffer),infile_ptr) != NULL){
		if(!strncmp(&buffer[0],"ATOM  ",6)){

			strncpy(res,&buffer[17],3);
			res[3]='\0';
            
            if(is_natural_amino(res)){
                n_aaa++;
            }else if(is_natural_nucleic(res)){
                n_naa++;
            }
            
            n_atoms++;
        }
    }
    
    CloseFile_B(&infile_ptr,"r");
    
    return n_naa / (float)n_atoms > 0.5f ? 1 : 0;
    
}