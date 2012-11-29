#include "flexaid.h"
#include "boinc.h"

#define NAA 20
#define NNA 4
#define NLIST 40

static char nucleic_acids[NNA][4]     = { "  A","  U","  C","  G" };
static char nucleic_acids_rev[NNA][4] = { "A  ","U  ","C  ","G  " };

static char protein_amino[NAA][4]     = { "GLY","ALA","VAL","LEU","ILE","MET",
					  "ASN","PRO","CYS","SER","THR","GLN",
					  "ASP","GLU","LYS","ARG","HIS",
					  "PHE","TRP","TYR" };


static char protein_atoms_order[NLIST][5]   = { " N  "," CA "," C  "," O  "," CB ",
						" CG "," SG "," OG "," CG1"," OG1"," CG2",
						" CD "," SD "," OD1"," CD1"," ND1"," OD2"," ND2"," CD2",
						" CE "," NE "," CE1"," NE1"," OE1"," OE2"," NE2"," CE2"," CE3",
						" CZ "," NZ "," CZ1"," CZ2"," CZ3",
						" CH "," CH1"," CH2"," OH "," NH1"," NH2",
						" OXT" };


void modify_pdb(char* infile, char* outfile, int exclude_het, int remove_water, int is_protein)
{
	char bufnul[10];
	char buffer[100];   // pdb line

	char lines[50][100]; // store residue lines
	int  nlines=0;

	int prev_resnum = -1;
	int resnum = -1;
	char res[4];

	int read = 0;
	int wrote = 0;

	FILE* infile_ptr = NULL;
	FILE* outfile_ptr = NULL;

	char insert = '-', prev_insert = '-';
	
	printf("Protein PDB files are reordered\n");
	printf("Hydrogens are removed\n");
	printf("'A' alternate conformation ONLY is chosen\n");
	printf("heterogroups are%s excluded\n", exclude_het ? "":" not");
	printf("water molecules will%s be removed\n", exclude_het || remove_water ? "":" not");
	
	if(!OpenFile_B(infile,"r",&infile_ptr)){
		fprintf(stderr,"ERROR: Could not order PDB file %s\n", infile);
		Terminate(20);
	}

	outfile_ptr = fopen(outfile,"w");
	if(outfile_ptr == NULL){
		fprintf(stderr, "ERROR: Could not write temporary PDB file.\n");
		Terminate(20);
	}
	
	
	while(fgets(buffer,sizeof(buffer),infile_ptr) != NULL){
		if(!strncmp(&buffer[0],"ATOM  ",6)){
			// all lines that start with 'ATOM  ' field

			read++;
			
			//0         1         2         3         4         5         6         
			//0123456789012345678901234567890123456789012345678901234567890123456789
			//ATOM     47  CB  ILE A   7      38.324  -3.725  17.587  1.00  0.00           C  
			strncpy(res,&buffer[17],3);
			res[3]='\0';

			strncpy(bufnul,&buffer[22],4);
			sscanf(bufnul,"%d",&resnum);
				
			// insertion of residue
			insert = buffer[26];
				
			// skip alternate conformations other than 'A'
			if(buffer[16] != ' ' && buffer[16] != 'A'){ continue; }
				
			if(resnum == prev_resnum && insert == prev_insert){
				if(is_protein && is_natural_amino(res)){
					// store line
					strcpy(lines[nlines++],buffer);
				}else if(!is_protein && is_natural_nucleic(res)){
					fprintf(outfile_ptr,"%s",buffer);
				}else{
					// ligands/mod. amino acids are marked as HETATM by default
					fprintf(outfile_ptr,"HETATM%s",&buffer[6]);
				}
					
			}else if(prev_resnum != -1){
				if(is_protein && nlines > 0){
					//write out ordered lines
					rewrite_residue2(lines,nlines,&wrote,outfile_ptr);
					nlines=0;
				}
				
				if(is_protein && is_natural_amino(res)){
					strcpy(lines[nlines++],buffer);
				}else if(!is_protein && is_natural_nucleic(res)){
					fprintf(outfile_ptr,"%s",buffer);					
				}else{
					// ligands/mod. amino acids are marked as HETATM by default
					fprintf(outfile_ptr,"HETATM%s",&buffer[6]);
				}
				
			}else{
				if(is_protein && is_natural_amino(res)){
					strcpy(lines[nlines++],buffer);
				}else if(is_natural_nucleic(res)){
					fprintf(outfile_ptr,"%s",buffer);
				}else{
					// ligands/mod. amino acids are marked as HETATM by default
					fprintf(outfile_ptr,"HETATM%s",&buffer[6]);
				}
			}
				
			prev_resnum = resnum;
			prev_insert = insert;
			

		}else{
			// all other lines that do not start with 'ATOM  ' field

			if(!strncmp(&buffer[0],"HETATM",6)){
				if(exclude_het) { continue; }
				else {
					if(!strncmp(&buffer[17],"HOH",3) && remove_water){ continue; }
				}
			}

			if(is_protein && nlines > 0){
				rewrite_residue2(lines,nlines,&wrote,outfile_ptr);
				nlines=0;
			}
			
			fprintf(outfile_ptr,"%s",buffer);
			
		}


	}

	if(is_protein && nlines > 0){	
		rewrite_residue2(lines,nlines,&wrote,outfile_ptr);
	}
	
	CloseFile_B(&infile_ptr,"r");

	fclose(outfile_ptr);
	
	printf("number of ATOM lines read is %d\n", read);
	printf("number of lines outputted is %d\n", wrote);
	
}

int get_NextLine(char lines[][100], int nlines){

	int k,l=-1,m=NLIST;
	char name[5];
	
	for(int i=0; i<nlines; i++){
		strncpy(name,&lines[i][12],4);
		name[4]='\0';
		
		k=NLIST;
		for(int j=NLIST-1; j>=0; --j)
			if(!strcmp(protein_atoms_order[j],name))
				k=j;
		
		if(k<m){
			m=k;
			l=i;
		}
			
	}

	return l;
}

void rewrite_residue2(char lines[][100], int nlines, int* wrote, FILE* outfile_ptr){
	
	int i;
	
	while((i=get_NextLine(lines,nlines)) != -1){
		fprintf(outfile_ptr,"%s",lines[i]);
		strcpy(lines[i],"                    ");
		(*wrote)++;
	}
	
}

int is_natural_amino(char* res){
	
	for(int i=0; i<NAA; i++){
		if(!strcmp(res,protein_amino[i])){
			return 1;
		}
	}

	return 0;
}

int is_natural_nucleic(char* res){
	
	for(int i=0; i<NNA; i++){
		if(!strcmp(res,nucleic_acids[i]) || !strcmp(res,nucleic_acids_rev[i])){
			return 1;
		}
	}

	return 0;
}
