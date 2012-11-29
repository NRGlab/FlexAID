#include "flexaid.h"
#include "boinc.h"

/*****************************************************************************
 * SUBROUTINE read_lig reads the ligand information from the input file 
 *****************************************************************************/
void read_lig(FA_Global* FA,atom** atoms,resid** residue,char ligfile[]){

	FILE *infile_ptr;        /* pointer to PDB file */
	char buffer[MAX_PATH__*2]; /* input line from PDB file */
	char field[7];           /* 6 letter field name, e.g. HETATM or ATOM */
	//char input_file[20];   /* input file name */
	char bufnul[30];         /* dumb string to read input strings */
	int  i,j,k,l;            /* dumb counter */
	int  flag;               /* a simple flag */
	char icfile[MAX_PATH__];   /* IC filename */
	//char anam[5];          /* temporary atom name */
	char rnam[4];            /* temporary residue name */
	//int  type;             /* temporary atom type */
	//char num_char[7];      /* string to read integer */
	//int  rec[3];           /* temporary array for atom reconstruction data */
	//int  opt[2];           /* temporary array for reading optimization data*/
	//char chain;
	int  list[MAX_ATM_HET];  /* array for atoms to be reconstructed */
	int  natm;               /* total number of entries in list */
	char help[10];
	int fdih;
	int nfdih;

	int bondlist[MAX_ATM_HET];
	int neighbours[MAX_ATM_HET];
	int nbonded;

	int **altfdih; //array used to set shift atoms
	char* pch;
	//int prot_latm_num;  
  
 
	FA->optres = (OptRes*)malloc(FA->MIN_OPTRES*sizeof(OptRes));
	if(!FA->optres){
		fprintf(stderr,"ERROR: Could not allocate memory for optres (ligand)\n");
		Terminate(2);
	}
	FA->MIN_OPTRES++;

	flag=0;
	FA->num_het=0;
	FA->num_het_atm=0;


	//prot_latm_num=(*atoms)[atm_cnt].number;
	//printf("prot_latm_num=%d\n",prot_latm_num);
	//PAUSE;
	//printf("ligfile=<%s>\n",ligfile);
	//PAUSE;

	infile_ptr=NULL;
	if (!OpenFile_B(ligfile,"r",&infile_ptr))
		Terminate(8);
 
  
	strcpy(icfile,ligfile);
	pch = strstr(icfile,".inp");
  
	if(pch != NULL){
		strcpy(pch,".ic");
		//printf("icfile = %s\n", icfile);
	}else{
		pch = icfile;
		pch[0] = '\0';
	}
  

	altfdih = (int**)malloc(FA->MIN_FLEX_BONDS*sizeof(int*));
	if(!altfdih){
		fprintf(stderr,"ERROR: memory allocation failed for altfdih.\n");
		Terminate(2);
	}
	for(i=0;i<FA->MIN_FLEX_BONDS;i++){
		altfdih[i] = (int*)malloc(3*sizeof(int));
		if(!altfdih[i]){
			fprintf(stderr,"ERROR: memory allocation failed for altfdih[%d].\n",i);
			Terminate(2);
		}
	}


	while (fgets(buffer, sizeof(buffer),infile_ptr)!=NULL){
    
		/*
		  01234567890123456789012345678901234567890123456789
		  HETTYP9000210  O12 m 900009000390004    0
		*/
		for (i=0;i<6;i++){field[i]=buffer[i];}
		field[6]='\0';
    
		// set residue.name residue.number and residue.chn
		if(strcmp(field,"RESIDU") == 0){
			for(i=0;i<=2;i++){rnam[i]=buffer[i+7];}
			rnam[3]='\0';
			FA->res_cnt++;
      
			if(FA->res_cnt==FA->MIN_NUM_RESIDUE){
				//printf("reallocating memory for residue\n");

				FA->MIN_NUM_RESIDUE++;

				(*residue) = (resid*)realloc((*residue),FA->MIN_NUM_RESIDUE*sizeof(resid));
				if(!(*residue)) {
					fprintf(stderr,"ERROR: memory allocation error for residue.\n");
					Terminate(2);
				}

				memset(&(*residue)[FA->MIN_NUM_RESIDUE-1],0,sizeof(residue));
				//printf("memory reallocated for residue\n");
			}
      
			(*residue)[FA->res_cnt].fatm = (int*)malloc(FA->MIN_ROTAMER*sizeof(int));
			(*residue)[FA->res_cnt].latm = (int*)malloc(FA->MIN_ROTAMER*sizeof(int));
			(*residue)[FA->res_cnt].bond = (int*)malloc(FA->MIN_FLEX_BONDS*sizeof(int));
			if (!(*residue)[FA->res_cnt].fatm ||
			    !(*residue)[FA->res_cnt].latm ||
			    !(*residue)[FA->res_cnt].bond){
				fprintf(stderr,"ERROR: memory allocation error for residue.fatm || .latm || .bond.\n");
				Terminate(2);
			}
      
			memset((*residue)[FA->res_cnt].fatm,0,FA->MIN_ROTAMER*sizeof(int));
			memset((*residue)[FA->res_cnt].latm,0,FA->MIN_ROTAMER*sizeof(int));
			memset((*residue)[FA->res_cnt].bond,0,FA->MIN_FLEX_BONDS*sizeof(int));
      
			FA->num_het++;
			FA->het_res[FA->num_het]=FA->res_cnt;
			(*residue)[FA->res_cnt].bonded=NULL;

			(*residue)[FA->res_cnt].type=1;
			strcpy((*residue)[FA->res_cnt].name,rnam);

			(*residue)[FA->res_cnt].chn=buffer[11];
			if((*residue)[FA->res_cnt].chn == '-'){(*residue)[FA->res_cnt].chn=' ';}

			(*residue)[FA->res_cnt].rot=0;
			(*residue)[FA->res_cnt].fdih=0;
			for(i=0;i<=3;i++){bufnul[i]=buffer[i+13];}
			bufnul[4]='\0';
			sscanf(bufnul,"%d",&(*residue)[FA->res_cnt].number);
		}
    
		//printf("residue[%d].rot = %d\n", res_cnt, (*residue)[res_cnt].rot);
    
		else if(strcmp(field,"HETTYP") == 0){
			FA->atm_cnt++;
			FA->atm_cnt_real++;
      
			if(FA->atm_cnt==FA->MIN_NUM_ATOM){
	
				//printf("reallocating memory for atoms\n");
	
				FA->MIN_NUM_ATOM += 50;
	
				(*atoms) = (atom*)realloc((*atoms),FA->MIN_NUM_ATOM*sizeof(atom));
				if(!(*atoms)){
					fprintf(stderr,"ERROR: memory allocation error for atoms.\n");
					Terminate(2);
				}

				memset(&(*atoms)[FA->MIN_NUM_ATOM-50],0,50*sizeof(atom));

				//printf("memory reallocated for atoms\n");
			}
      

			//////////////////////////////////
			// num_atm assignment

			for (i=0;i<5;i++){help[i]=buffer[6+i];}
			help[5]='\0';
			sscanf(help,"%d",&i);
			(*atoms)[FA->atm_cnt].number=i;
      
			FA->num_atm[i]=FA->atm_cnt;

			////////////////////////////////////

      
			(*atoms)[FA->atm_cnt].bond[0]=0;
			//printf("atm_cnt=%d\n",FA->atm_cnt);
			// set reside.fatm and residue.latm
			if(flag == 0){
				(*residue)[FA->res_cnt].fatm[0]=FA->atm_cnt;

				/*
				  printf("FA->het_res[%d]=%d\t(*residue)[FA->het_res[%d]].fatm[0]=%d\n",
				  FA->num_het,
				  FA->het_res[FA->num_het],
				  FA->num_het,
				  (*residue)[FA->het_res[FA->num_het]].fatm[0]);
				*/

				flag=1;
			}
			(*residue)[FA->res_cnt].latm[0]=FA->atm_cnt;
      
			//HETTYP _____ _ ____ _ __________________
      
			help[0]=buffer[11];
			help[1]=buffer[12];
			help[2]='\0';
			sscanf(help,"%d",&(*atoms)[FA->atm_cnt].type);
      
			for (i=0;i<4;i++){(*atoms)[FA->atm_cnt].name[i]=buffer[14+i];}
			(*atoms)[FA->atm_cnt].name[4]='\0';
      
			if((*atoms)[FA->atm_cnt].type > FA->ntypes){
				printf("WARNING: res %s atom %s has atom type %d when %d types are defined\n",
				       (*residue)[FA->res_cnt].name, (*atoms)[FA->atm_cnt].name, (*atoms)[FA->atm_cnt].type, FA->ntypes);
				printf("WARNING: type %d is set to neutral (6)\n", (*atoms)[FA->atm_cnt].type);
	      
				(*atoms)[FA->atm_cnt].type = 6;
			}

			(*atoms)[FA->atm_cnt].recs=buffer[19];
			(*atoms)[FA->atm_cnt].ncons=0;
			(*atoms)[FA->atm_cnt].cons=NULL;
			(*atoms)[FA->atm_cnt].graph=0;
			(*atoms)[FA->atm_cnt].coor_ref=NULL;
			
			for (j=0;j<3;j++){
				for (i=0;i<5;i++){
					help[i]=buffer[21+i+j*5];
				}
				help[5]='\0';
				sscanf(help,"%d",&(*atoms)[FA->atm_cnt].rec[j]);
			}
      
			for(i=0;i<5;i++){
				help[i]=buffer[36+i];
			}
			help[5]='\0';
			sscanf(help,"%d",&(*atoms)[FA->atm_cnt].graph);


			// set num_atm[] for the ligand
			FA->num_atm[(*atoms)[FA->atm_cnt].number]=FA->atm_cnt;
      
			// set atoms.radius
			(*atoms)[FA->atm_cnt].radius=assign_radius((*atoms)[FA->atm_cnt].name);
      
			// set atoms.ofres
			(*atoms)[FA->atm_cnt].ofres=FA->res_cnt;
      
			/* printf("(%s) (%d) (%d) (%s) (%c) (%d) (%d) (%d) (%d)\n",
			   field,
			   (*atoms)[FA->atm_cnt].number,
			   (*atoms)[FA->atm_cnt].type,
			   (*atoms)[FA->atm_cnt].name,
			   (*atoms)[FA->atm_cnt].recs,
			   (*atoms)[FA->atm_cnt].rec[0],
			   (*atoms)[FA->atm_cnt].rec[1],
			   (*atoms)[FA->atm_cnt].rec[2],
			   (*atoms)[FA->atm_cnt].rec[3]
			   );
			*/
			FA->num_het_atm++;
		}
    
		// set residue.fdih and residue.bond
		else if(strcmp(field,"FLEDIH") == 0){ 
			(*residue)[FA->res_cnt].fdih++;
			//printf("fdih=%d\n",(*residue)[FA->res_cnt].fdih);

			if((*residue)[FA->res_cnt].fdih==FA->MIN_FLEX_BONDS){
				FA->MIN_FLEX_BONDS += 5;
	      
				// fdih
				//printf("re-allocating memory for fdih...\n");
				(*residue)[FA->res_cnt].bond = (int*)realloc((*residue)[FA->res_cnt].bond,FA->MIN_FLEX_BONDS*sizeof(int));
				if(!(*residue)[FA->res_cnt].bond){
					fprintf(stderr,"ERROR: memory allocation error for residue.bond.\n");
					Terminate(2);
				}
				//printf("memory re-allocated fdih...\n");
				memset(&(*residue)[FA->res_cnt].bond[FA->MIN_FLEX_BONDS-5],0,5*sizeof(int));
	      
	      
				// altfdih
				altfdih = (int**)realloc(altfdih,FA->MIN_FLEX_BONDS*sizeof(int*));
				if(!altfdih){
					fprintf(stderr,"ERROR: memory allocation failed for altfdih.\n");
					Terminate(2);
				}
	      
				for(i=FA->MIN_FLEX_BONDS-5;i<FA->MIN_FLEX_BONDS;i++){
					altfdih[i] = (int*)malloc(3*sizeof(int));
					if(!altfdih[i]){
						fprintf(stderr,"ERROR: memory allocation failed for altfdih[%d].\n",i);
						Terminate(2);
					}
				}
			}
      
			for(k=0;k<3;k++) altfdih[(*residue)[FA->res_cnt].fdih][k]=0;
      
			for(i=7;i<9;i++) bufnul[i-7]=buffer[i];
			bufnul[2]='\0';
			sscanf(bufnul,"%d",&fdih);
      
			nfdih=(int)(strlen(buffer)-11)/5;
			//printf("nfdih=%d\n",nfdih);
			for(i=0;i<nfdih;i++){
				for(j=0;j<5;j++){
					bufnul[j]=buffer[i*5+j+10];
				}
				bufnul[5]='\0';
				sscanf(bufnul,"%d",&k);
	
				altfdih[(*residue)[FA->res_cnt].fdih][i]=FA->num_atm[k];

				if ((*residue)[FA->res_cnt].bond[(*residue)[FA->res_cnt].fdih]==0)
					(*residue)[FA->res_cnt].bond[(*residue)[FA->res_cnt].fdih]=FA->num_atm[k];

			        //printf("added %d to bond=%d\n",(*residue)[FA->res_cnt].bond[(*residue)[FA->res_cnt].fdih],(*residue)[FA->res_cnt].fdih);
			}
      
			//old code
			//sscanf(buffer,"%s %d %d",field,&i,&l);
			//(*residue)[FA->res_cnt].bond[(*residue)[FA->res_cnt].fdih]=FA->num_atm[l];
		        //printf("res_cnt: %d\tbond #(fdih): %d\tnum_atm: %d\n", FA->res_cnt, (*residue)[FA->res_cnt].fdih, FA->num_atm[l]);
		}
    
		// set residue.gpa
		else if(strcmp(field,"GPATOM") == 0){
			// only allocate memory for gpa for ligand
			(*residue)[FA->res_cnt].gpa = (int*)malloc(3*sizeof(int));
			if(!(*residue)[FA->res_cnt].gpa){
				fprintf(stderr,"ERROR: memory allocation error for residue.gpa.\n");
				Terminate(2);
			}
      
			sscanf(buffer,"%s %d %d %d",field,&i,&j,&l);
			/* printf("i: %d j: %d l: %d\n",i,j,l);
			   printf("num_atm[i]: %d num_atm[j]: %d num_atm[l]: %d\n",
			   FA->num_atm[i],
			   FA->num_atm[j],
			   FA->num_atm[l]);*/
			(*residue)[FA->res_cnt].gpa[0]=FA->num_atm[i];
			(*residue)[FA->res_cnt].gpa[1]=FA->num_atm[j];
			(*residue)[FA->res_cnt].gpa[2]=FA->num_atm[l];
			//printf("res_cnt[%d]=gpa0: %d gpa1: %d gpa2: %d\n",FA->res_cnt,(*residue)[FA->res_cnt].gpa[0],(*residue)[FA->res_cnt].gpa[1],(*residue)[FA->res_cnt].gpa[2]);
		}

    
		// set atoms.bond
		else if(strcmp(field,"CONECT") == 0){
			read_conect(FA,atoms,buffer);
		}
    
		/*
		// reads IC data file name
		if(strcmp(field,"ICDATA") == 0){
		sscanf(buffer,"%s %s",field,icfile);
		//printf("BUFFER=<%s>\n",buffer);
		//printf("ICFILE=<%s>\n",icfile);
		//PAUSE;
		}
		*/

		else if(strcmp(field,"ENDINP") == 0){ 
			CloseFile_B(&infile_ptr,"r");
			break;
		}
	}
  
  
	// corrects the value of rec[] to internal number scheme
	for(i=(*residue)[FA->res_cnt].fatm[0];i<=(*residue)[FA->res_cnt].latm[0];i++){
		if((*atoms)[i].recs == 'm'){
      
			for(l=0;l<=3;l++){
				if((*atoms)[i].rec[l] != 0){
					//printf("num_atm[%d]=%d ",(*atoms)[i].rec[l],FA->num_atm[(*atoms)[i].rec[l]]);
					(*atoms)[i].rec[l]=FA->num_atm[(*atoms)[i].rec[l]];
				}
	
				//printf("atoms[%d].rec[%d]=%d ",FA->atm_cnt,l,(*atoms)[FA->atm_cnt].rec[l]);
			}
			//PAUSE;
		}
	}
  
	infile_ptr=NULL;
	// reads icfile to set atoms.dis atoms.ang and atoms.dih
	if (!OpenFile_B(icfile,"r",&infile_ptr))
		Terminate(8);
  
	while (fgets(buffer, sizeof(buffer),infile_ptr)){

		for (i=0;i<6;i++){field[i]=buffer[i];}
		field[6]='\0';
	  
		// set residue.reference pcg
		if(strcmp(field,"REFPCG") == 0){      
			for(i=0;i<=25;i++){bufnul[i]=buffer[7+i];}
			bufnul[26]='\0';
			sscanf(bufnul,"%f %f %f",&FA->ori[0],&FA->ori[1],&FA->ori[2]);

		}else{

			for(i=0;i<=4;i++){bufnul[i]=buffer[i];}
			bufnul[5]='\0';
			sscanf(bufnul,"%d",&l);
			l=FA->num_atm[l];
			//printf("l=<%d>\n",l);
			for(i=0;i<=25;i++){bufnul[i]=buffer[7+i];}
			bufnul[26]='\0';
			sscanf(bufnul,"%f %f %f",&(*atoms)[l].dis,&(*atoms)[l].ang,&(*atoms)[l].dih);
		  
			(*atoms)[l].dis += 0.0001f;
			(*atoms)[l].ang += 0.0001f;
			(*atoms)[l].dih += 0.0001f;
			//printf("%8.3f %8.3f %8.3f\n",(*atoms)[l].dis,(*atoms)[l].ang,(*atoms)[l].dih);
		}
	}
	CloseFile_B(&infile_ptr,"r");
	//PAUSE;

	printf("the protein center of coordinates is: %8.3f %8.3f %8.3f\n", FA->ori[0],FA->ori[1],FA->ori[2]);

	// creates list of atoms in ligand to create cc
	buildlist(FA,*atoms,*residue,FA->res_cnt,0,&natm,list);
	//  for(i=1;i<=residue[res_cnt].fdih;i++)
	//printf("%2d=%6d%6d%6d\n",i,altfdih[i][0],altfdih[i][1],altfdih[i][2]);

	if((*residue)[FA->res_cnt].fdih > 0)
		assign_shift(*atoms,*residue,FA->res_cnt,natm,list,altfdih);


	FA->optres[0].rnum=FA->res_cnt;
	FA->optres[0].type=1;
	FA->optres[0].tot=FA->num_het_atm;

	//*****NEW*****
	//neighbours n-bonds-away are stored in memory
	for(i=0;i<natm;i++){
		nbonded=0;
    
		bondedlist(*atoms,list[i],FA->bloops,&nbonded,bondlist,neighbours);
		update_bonded(&(*residue)[FA->res_cnt],natm,nbonded,bondlist,neighbours);

		/*
		  printf("bondlist[%d] for atom %d\t",i,(*atoms)[list[i]].number);
		  for(j=0;j<nbonded;j++) printf("%6d",(*atoms)[bondlist[j]].number);
		  printf("\n");
		*/
	}
  
	// prints bonded matrix
	/*
	  for(i=0;i<natm;i++){
	  printf("\t%5d\t",i);
	  for(j=0;j<natm;j++){
	  printf("%2d",(*residue)[FA->res_cnt].bonded[i][j]);
	  }
	  printf("\n");
	  }
	  getchar();
	*/

	FA->num_optres++;

	buildcc(FA,*atoms,natm,list);
	
	for(i=(*residue)[FA->res_cnt].fatm[0];i<=(*residue)[FA->res_cnt].latm[0];i++){
		//printf("%d %8.3f %8.3f %8.3f\n",(*atoms)[i].number,(*atoms)[i].coor[0],(*atoms)[i].coor[1],(*atoms)[i].coor[2]);
		for(j=0; j<3; j++){
			(*atoms)[i].coor_ori[j] = (*atoms)[i].coor[j];
		}
	}

	for(i=0; i<FA->MIN_FLEX_BONDS; i++){
		free(altfdih[i]);
	}
	free(altfdih);

	return;
}
