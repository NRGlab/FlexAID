#include "flexaid.h"
#include "boinc.h"

void read_emat(FA_Global* FA, char* scr_bin_file){
	
	int i,j,l;
	FILE* infile_ptr;
	float z;
	char buffer[40];
	char bufnul[15];

	infile_ptr=NULL;
	if (!OpenFile_B(scr_bin_file,"r",&infile_ptr)){
		fprintf(stderr,"ERROR: Could not read input file: %s\n", scr_bin_file);
		Terminate(8);
	}
	
	l=0;
	while(fgets(buffer, sizeof(buffer), infile_ptr) != NULL){ l++; }
	//if(strchr(buffer, '=') != NULL){ l++; }	  
	
	z = zero(1.0/2.0f, 1.0/2.0f, (float)(-l));
	if(fabs(z - (float)((int)z)) > 0.001) {
		fprintf(stderr,"ERROR: Number of items in %s is incorrect\n", scr_bin_file);
		Terminate(12);
	}
	
	FA->ntypes = (int)(z + 0.001);
	printf("number of atom types: %d\n", FA->ntypes);
	
	FA->energy = (float**)malloc((FA->ntypes+1)*sizeof(float*));
	if(FA->energy == NULL){
		fprintf(stderr,"ERROR: could not allocate memory for FA->energy\n");
		Terminate(2);
	}
	
	for(i=1; i<=FA->ntypes; i++){
		FA->energy[i] = (float*)malloc((FA->ntypes+1)*sizeof(float));
		if(FA->energy[i] == NULL){
			fprintf(stderr,"ERROR: could not allocate memory for FA->energy[%d]\n", i);
			Terminate(2);
		}
	}
	
	fseek(infile_ptr, 0, SEEK_SET);
	for(i=1;i<=FA->ntypes;i++){
		for(j=i;j<=FA->ntypes;j++){
			if(fgets(buffer, sizeof(buffer), infile_ptr) != NULL){
				//for(l=0;l<=7;l++){bufnul[l]=buffer[l+13];}
				//bufnul[8]='\0';
				//sscanf(bufnul,"%f",&FA->energy[i][j]);
				FA->energy[i][j] = atof(&buffer[13]);
				FA->energy[j][i]=FA->energy[i][j];
			}
		}
	}
	
	CloseFile_B(&infile_ptr,"r");
	
}
