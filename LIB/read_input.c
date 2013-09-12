#include "flexaid.h"
#include "boinc.h"

/***************************************************************************** 
 * SUBROUTINE read_input reads input file.
 *****************************************************************************/
void read_input(FA_Global* FA,atom** atoms, resid** residue,rot** rotamer,gridpoint** cleftgrid,char* input_file){

	FILE *infile_ptr;          /* pointer to input file */
	char buffer[MAX_PATH__*2];   /* a line from the INPUT file */
	//char input_file[20];     /* input file name */
	char field[7];             /* field names on INPUT file */
	int  i,k;              /* dumb counter */

	//int  flag;               /* a simple flag */
	char pdb_name[MAX_PATH__];       /* 4 letter PDB filename */
	char lig_file[MAX_PATH__];       /* 4 letter PDB filename */
	char rmsd_file[MAX_PATH__];     /* 4 letter PDB filename */
	char clf_file[MAX_PATH__];      /* Cleft file that contains sphere coordinates */
	char normal_file[MAX_PATH__];   /* normal mode grid file */
	char eigen_file[MAX_PATH__];   /* normal mode grid file */

	char emat_forced[MAX_PATH__];
	char emat[MAX_PATH__];      /* interaction matrix file*/

	char deftyp_forced[MAX_PATH__];
	char deftyp[MAX_PATH__];     /* amino/nucleotide definition file*/

	char constraint_file[MAX_PATH__]; /* constraint file */

	char rotlib_file[MAX_PATH__];  /* rotamer library file */
	char rotobs_file[MAX_PATH__];  /* rotamer observations file */

	char rngopt[7];
	char rngoptline[MAX_PATH__];

	//char anam[5];            /* temporary atom name */
	//char rnam[4];            /* temporary residue name */
	//int  type;               /* temporary atom type */
	//char num_char[7];        /* string to read integer */
	//int  rec[3];             /* temporary array for atom reconstruction data */
	int opt[2];              /* temporary array for reading optimization data*/
	char chain='-';
	char a[7],b[7]; //,mol_name[4]; 

	char flexscfile[MAX_PATH__];
	char gridfile[MAX_PATH__];
	char gridfilename[MAX_PATH__];
	//char sphfile[MAX_PATH__];
	char tmpprotname[MAX_PATH__];

	char optline[MAX_PAR][MAX_PATH__];

	// opt lines counter
	int nopt=0;

	sphere *spheres, * _sphere;
	
	flexscfile[0]='\0';
	rmsd_file[0]='\0';
	FA->state_path[0]='\0';
	FA->temp_path[0]='\0';
	constraint_file[0] = '\0';

	deftyp_forced[0] = '\0';
	emat_forced[0] = '\0';
	FA->dependencies_path[0] = '\0';
  
	spheres = NULL;

	infile_ptr=NULL;
	if (!OpenFile_B(input_file,"r",&infile_ptr)){
		fprintf(stderr,"ERROR: Could not read input file: %s\n",input_file);
		Terminate(8);
	}

	while (fgets(buffer, sizeof(buffer),infile_ptr)){

		buffer[strlen(buffer)-1] = '\0';
        
		for(i=0;i<6;++i) field[i]=buffer[i];
		field[6]='\0';
		
		if(strcmp(field,"PDBNAM") == 0){strcpy(pdb_name,&buffer[7]);}
		if(strcmp(field,"INPLIG") == 0){strcpy(lig_file,&buffer[7]);}
		if(strcmp(field,"METOPT") == 0){sscanf(buffer,"%s %s",a,FA->metopt);}
		if(strcmp(field,"DEEFLX") == 0){FA->deelig_flex=1;}
		if(strcmp(field,"BPKENM") == 0){sscanf(buffer,"%s %s",a,FA->bpkenm);}
		if(strcmp(field,"COMPLF") == 0){sscanf(buffer,"%s %s",a,FA->complf);}
		if(strcmp(field,"VCTSCO") == 0){sscanf(buffer,"%s %s",a,FA->vcontacts_self_consistency);}
		if(strcmp(field,"RNGOPT") == 0){strcpy(rngoptline,buffer);for(i=0;i<6;i++)rngopt[i]=buffer[7+i];rngopt[6]='\0';}
		if(strcmp(field,"OPTIMZ") == 0){strcpy(optline[nopt++],buffer);}
		if(strcmp(field,"DEFTYP") == 0){strcpy(deftyp_forced,&buffer[7]);}
		if(strcmp(field,"FLEXSC") == 0){strcpy(flexscfile,&buffer[7]);}
		if(strcmp(field,"NMAMOD") == 0){sscanf(buffer,"%s %d",a,&FA->normal_modes);}
		if(strcmp(field,"NMAAMP") == 0){strcpy(normal_file,&buffer[7]);}
		if(strcmp(field,"NMAEIG") == 0){strcpy(eigen_file,&buffer[7]);}
		if(strcmp(field,"RMSDST") == 0){strcpy(rmsd_file,&buffer[7]);}
		if(strcmp(field,"CLRMSD") == 0){sscanf(buffer,"%s %f",a,&FA->cluster_rmsd);}
		if(strcmp(field,"EXCHET") == 0){FA->exclude_het=1;}
		if(strcmp(field,"INCHOH") == 0){FA->remove_water=0;}
		if(strcmp(field,"PERMEA") == 0){sscanf(buffer,"%s %f",field,&FA->permeability);}
		if(strcmp(field,"INTRAF") == 0){sscanf(buffer,"%s %f",field,&FA->intrafraction);}
		if(strcmp(field,"VARDIS") == 0){sscanf(buffer,"%s %lf",field,&FA->delta_angstron);}
		if(strcmp(field,"VARANG") == 0){sscanf(buffer,"%s %lf",field,&FA->delta_angle);}
		if(strcmp(field,"VARDIH") == 0){sscanf(buffer,"%s %lf",field,&FA->delta_dihedral);}
		if(strcmp(field,"VARFLX") == 0){sscanf(buffer,"%s %lf",field,&FA->delta_flexible);}
		if(strcmp(field,"SLVPEN") == 0){sscanf(buffer,"%s %f",field,&FA->solventterm);}
		if(strcmp(field,"OUTRNG") == 0){FA->output_range=1;}
		if(strcmp(field,"USEDEE") == 0){FA->useflexdee=1;}
		if(strcmp(field,"IMATRX") == 0){strcpy(emat_forced,&buffer[7]);}
		if(strcmp(field,"DEECLA") == 0){sscanf(buffer,"%s %f",field,&FA->dee_clash);}
		if(strcmp(field,"CONSTR") == 0){strcpy(constraint_file,&buffer[7]);}
		if(strcmp(field,"MAXRES") == 0){sscanf(buffer,"%s %d",field,&FA->max_results);}
		if(strcmp(field,"SPACER") == 0){sscanf(buffer,"%s %f",field,&FA->spacer_length);}
		if(strcmp(field,"DEPSPA") == 0){strcpy(FA->dependencies_path,&buffer[7]);}
		if(strcmp(field,"STATEP") == 0){strcpy(FA->state_path,&buffer[7]);}
		if(strcmp(field,"TEMPOP") == 0){strcpy(FA->temp_path,&buffer[7]);}
		if(strcmp(field,"NRGSUI") == 0){FA->nrg_suite=1;}
		if(strcmp(field,"NRGOUT") == 0){sscanf(buffer,"%s %d",field,&FA->nrg_suite_timeout);}
		if(strcmp(field,"SCOOUT") == 0){FA->output_scored_only=1;}
	}
    
	CloseFile_B(&infile_ptr,"r");
  

	//////////////////////////////////////////////////
	////////// read input files afterwards ///////////
	//////////////////////////////////////////////////
  
	printf("dependencies path=%s\n", FA->dependencies_path);

	// default state path (controls pause-stop-abort)
	if(!strcmp(FA->state_path,"")){
		strcpy(FA->state_path,FA->base_path);
	}

	if(!strcmp(FA->temp_path,"")){
		strcpy(FA->temp_path,FA->base_path);
	}
    
	// temporary pdb name
	strcpy(tmpprotname,FA->temp_path);
    
#ifdef _WIN32
	strcat(tmpprotname,"\\target.pdb");
#else
	strcat(tmpprotname,"/target.pdb");
#endif


	/************************************************************/
	/********          INTERACTION MATRIX              **********/
	/************************************************************/
	if(!strcmp(emat_forced,"")){
		// use default
		if(!strcmp(FA->dependencies_path,"")){
			// use executable path
			strcpy(emat,FA->base_path);		
		}else{
			// use new dependencies path
			strcpy(emat,FA->dependencies_path);
		}
#ifdef _WIN32
		strcat(emat,"\\M6_cons_3.dat");
#else
		strcat(emat,"/M6_cons_3.dat");
#endif
	}else{
		// use forced matrix
		strcpy(emat,emat_forced);
	}
	
	printf("interaction matrix is <%s>\n", emat);
	read_emat(FA,emat);

	printf("pdb target is <%s>\n", pdb_name);
	if(rna_structure(pdb_name)){
		printf("target molecule is a RNA structure\n");
		FA->is_protein = 0;
	}
	
	/************************************************************/
	/********          DEFINITION OF TYPES             **********/
	/************************************************************/
	if(!strcmp(deftyp_forced,"")){
		// use default definition
		
		if(!strcmp(FA->dependencies_path,"")){
			strcpy(deftyp,FA->base_path);
		}else{
			strcpy(deftyp,FA->dependencies_path);
		}
		
		
		if(FA->is_protein){
#ifdef _WIN32
			strcat(deftyp,"\\AMINO.def");
#else
			strcat(deftyp,"/AMINO.def");
#endif
		}else{
#ifdef _WIN32
			strcat(deftyp,"\\NUCLEOTIDES.def");
#else
			strcat(deftyp,"/NUCLEOTIDES.def");
#endif
		}
		
	}else{
		// use forced definition of types
		strcpy(deftyp,deftyp_forced);
	}

	printf("definition of types is <%s>\n", deftyp);

	///////////////////////////////////////////////
  
	printf("read PDB file <%s>\n",pdb_name);

	modify_pdb(pdb_name,tmpprotname,FA->exclude_het,FA->remove_water,FA->is_protein);
	read_pdb(FA,atoms,residue,tmpprotname);
	
	(*residue)[FA->res_cnt].latm[0]=FA->atm_cnt;
	for(k=1;k<=FA->res_cnt;k++){
		FA->atm_cnt_real += (*residue)[k].latm[0]-(*residue)[k].fatm[0]+1;
	}
  
	calc_center(FA,*atoms,*residue);
  
	if(FA->is_protein){ residue_conect(FA,*atoms,*residue,deftyp); }

	//assign_types(FA,*atoms,*residue,deftyp);

	//////////////////////////////////////////////

	printf("read ligand input file <%s>\n",lig_file);
	read_lig(FA,atoms,residue,lig_file);

	//////////////////////////////////////////////
	
	assign_radii_types(FA,(*atoms),(*residue));    
	printf("radii are now assigned\n");
    
	if(strcmp(rmsd_file,"")){
		printf("read rmsd structure <%s>: will match atom numbers\n",rmsd_file);

		int rmsd_atoms = read_rmsdst(FA,*atoms,*residue,rmsd_file);

		if(rmsd_atoms){
			FA->refstructure=1;
			printf("will use %d atoms to calculate RMSD\n", rmsd_atoms);
		}else{
			printf("no atoms is used to calculate RMSD\n");
		}
	}

	//////////////////////////////////////////////
  
	if(FA->normal_modes > 0){
		printf("read files related to NMA\n");
		read_normalgrid(FA,normal_file);
		read_eigen(FA,eigen_file);
		assign_eigen(FA,*atoms,*residue,FA->res_cnt,FA->normal_modes);
	}

  
	//////////////////////////////////////////////

	if(strcmp(flexscfile,"") && FA->is_protein){
		
		printf("read flexsc file <%s>\n",flexscfile);
		read_flexscfile(FA,*residue,rotamer,flexscfile,rotlib_file,rotobs_file);

		
		if (FA->rotobs) {

			// use rotamers found in bound conformations (hap2db)
			if(!strcmp(FA->dependencies_path,"")){
				strcpy(rotobs_file,FA->base_path);
			}else{
				strcpy(rotobs_file,FA->dependencies_path);
			}
            
#ifdef _WIN32
			strcat(rotobs_file,"\\rotobs.lst");
#else
			strcat(rotobs_file,"/rotobs.lst");
#endif
            
			printf("read rotamer observations <%s>\n",rotobs_file);
			read_rotobs(FA,rotamer,rotobs_file);
		}else{
			
			// use penultimate rotamer library instances
			if(!strcmp(FA->dependencies_path,"")){
				strcpy(rotlib_file,FA->base_path);
			}else{
				strcpy(rotlib_file,FA->dependencies_path);
			}
            
#ifdef _WIN32
			strcat(rotlib_file,"\\Lovell_LIB.dat");
#else
			strcat(rotlib_file,"/Lovell_LIB.dat");
#endif
            
			printf("read rotamer library <%s>\n",rotlib_file);
			read_rotlib(FA,rotamer,rotlib_file);
		}
    
		if(FA->nflxsc > 0 && FA->rotlibsize > 0){
			build_rotamers(FA,atoms,*residue,*rotamer); 
			//build_close(FA,residue,atoms);          
		}

	}

	//////////////////////////////////////////////

	if(strcmp(constraint_file,"")){
		printf("read constraint_file <%s>\n", constraint_file);
		read_constraints(FA,*atoms,*residue,constraint_file);
		assign_constraint_threshold(FA,*atoms,FA->constraints,FA->num_constraints);

	}

	///////////////////////////////////////////////

	if(!strcmp(rngopt,"LOCCEN")){
		strcpy(FA->rngopt,"loccen");
		
		_sphere = (sphere*)malloc(sizeof(sphere));
		if(_sphere == NULL){
			fprintf(stderr,"ERROR: memory allocation error for spheres (LOCCEN)\n");
			Terminate(2);
		}
        
		sscanf(rngoptline,"%s %s %f %f %f %f",
		       a,b,
		       &_sphere->center[0],&_sphere->center[1],&_sphere->center[2],
		       &_sphere->radius);
		_sphere->prev = NULL;
		
		// point to the new sphere created
		spheres = _sphere;
		
		(*cleftgrid) = generate_grid(FA,spheres,(*atoms),(*residue));
		calc_cleftic(FA,*cleftgrid);
        
	}else if(!strcmp(rngopt,"LOCCLF")){
        
		//RNGOPT LOCCLF filename.pdb
		strcpy(FA->rngopt,"locclf");
		strcpy(clf_file,&rngoptline[14]);
        
		printf("read binding-site grid <%s>\n",clf_file);
		spheres = read_spheres(clf_file);
        
		(*cleftgrid) = generate_grid(FA,spheres,(*atoms),(*residue));
		calc_cleftic(FA,*cleftgrid);
	}
    
	//printf("IC bounds...\n");
	ic_bounds(FA,FA->rngopt);

	///////////////////////////////////////////////////////////////////////////////
    
	for(i=0;i<nopt;i++){
		//printf("%s\n", optline[i]);
		if(i==MAX_PAR){
			printf("WARNING: number of params allowed was reached (100). other params will be skipped.\n");
			break;
		}
    
		//printf("optline[%d]=%s\n",i,optline[i]);

		sscanf(optline[i],"%s %d %s %d",a,&opt[0],a,&opt[1]);
		//printf("%d %d\n",opt[0],opt[1]);
		//getchar();
		//chain=buffer[11];
		chain=a[0];
		if(chain == '-'){chain = ' ';}
		//printf("Add2 optimiz vector...\n");
		add2_optimiz_vec(FA,*atoms,*residue,opt,chain,"");

	}

	add2_optimiz_vec(FA,*atoms,*residue,opt,chain,"SC");
	add2_optimiz_vec(FA,*atoms,*residue,opt,chain,"NM");
    
	//////////////////////////////////////////////
    
	if(FA->translational && strcmp(rngopt,"LOCCEN") && strcmp(rngopt,"LOCCLF")){
		fprintf(stderr,"ERROR: the binding-site is not defined\n");
		Terminate(2);
	}
    
	if(FA->output_range){		
        
#ifdef _WIN32
		strcpy(gridfilename,"\\grid.sta.pdb");
#else
		strcpy(gridfilename,"/grid.sta.pdb");
#endif
        
		strcpy(gridfile,FA->temp_path);
		strcat(gridfile,gridfilename);
        
		write_grid(FA,*cleftgrid,gridfile);
	}
    
	if(FA->translational && FA->num_grd==1){
		fprintf(stderr,"ERROR: the binding-site has no anchor points\n");
		Terminate(2);
	}
    
	// fill in optres pointer in atoms struct.
	update_optres(*atoms,FA->atm_cnt,FA->optres,FA->num_optres);
  	
	if(FA->nrg_suite){
		if(FA->translational){
			for(i=1; i<FA->num_grd; i++){
				printf("Grid[%d]=%8.3f%8.3f%8.3f\n", i, (*cleftgrid)[i].coor[0], (*cleftgrid)[i].coor[1], (*cleftgrid)[i].coor[2]);
				if(i % 1000 == 0){
					fflush(stdout);
				}
				fflush(stdout);
			}
		}
	}
    
	/* FREE SPHERES LINKED-LIST */
	while(spheres != NULL){
		_sphere = spheres->prev;
		free(spheres);

		spheres = _sphere;
	}

	return;
}
