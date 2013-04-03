#include "gaboom.h"
#include "boinc.h"
#include "Vcontacts.h"

int main(int argc, char **argv){

	int   i,j;
	int   natm;

	char remark[MAX_REMARK];
	char tmpremark[MAX_REMARK];
	char dockinp[MAX_PATH__];
	char gainp[MAX_PATH__];
	char *pch;                               // for finding base path
	char end_strfile[MAX_PATH__];
	char tmp_end_strfile[MAX_PATH__];

	int memchrom=0;
  
	time_t sta_timer,end_timer;
	struct tm *sta,*end;
	int sta_val[3],end_val[3];
	long ct; // computational time

	atom *atoms = NULL;
	resid *residue = NULL;
	resid *res_ptr = NULL;
	cfstr cf;
	cfstr* cf_ptr = NULL;
	rot* rotamer = NULL;
	chromosome* chrom = NULL;
	chromosome* chrom_snapshot = NULL;
	genlim* gene_lim = NULL;
	gridpoint* cleftgrid = NULL;

	//flexaid global variables
	FA_Global* FA = NULL;
	GB_Global* GB = NULL;
	VC_Global* VC = NULL;

	FA = (FA_Global*)malloc(sizeof(FA_Global));
	GB = (GB_Global*)malloc(sizeof(GB_Global));
	VC = (VC_Global*)malloc(sizeof(VC_Global));

	if(!FA || !GB || !VC){
		fprintf(stderr,"ERROR: Could not allocate memory for FA || GB || VC\n");
		Terminate(2);
	}

	memset(FA,0,sizeof(FA_Global));
	memset(GB,0,sizeof(GB_Global));
	memset(VC,0,sizeof(VC_Global));

	FA->contacts = (int*)malloc(100000*sizeof(int));
	if(FA->contacts == NULL){
		fprintf(stderr,"ERROR: Could not allocate memory for contacts\n");
		Terminate(2);
	}

	VC->ptorder = (ptindex*)malloc(MAX_PT*sizeof(ptindex));
	VC->centerpt = (vertex*)malloc(MAX_PT*sizeof(vertex));
	VC->poly = (vertex*)malloc(MAX_POLY*sizeof(vertex));
	VC->cont = (plane*)malloc(MAX_PT*sizeof(plane));
	VC->vedge = (edgevector*)malloc(MAX_POLY*sizeof(edgevector));

	if(!VC->ptorder || !VC->centerpt || !VC->poly ||
	   !VC->cont || !VC->vedge){
		fprintf(stderr,"ERROR: Could not allocate memory for ptorder || centerpt || poly || cont || vedge\n");
		Terminate(2);
	}


	// set minimal default values
	FA->MIN_NUM_ATOM = 1000;
	FA->MIN_NUM_RESIDUE = 250;
	FA->MIN_ROTAMER_LIBRARY_SIZE = 155;
	FA->MIN_ROTAMER = 1;
	FA->MIN_FLEX_BONDS = 5;
	FA->MIN_CLEFTGRID_POINTS = 250;  
	FA->MIN_PAR = 6;  
	FA->MIN_FLEX_RESIDUE = 5;  
	FA->MIN_NORMAL_GRID_POINTS = 250;
	FA->MIN_OPTRES = 1;
	FA->MIN_CONSTRAINTS = 1;

	FA->num_optres = 0;

	FA->normal_grid = NULL;
	FA->supernode = 0;
	FA->eigenvector = NULL;
	FA->psFlexDEENode = NULL;
	FA->FlexDEE_Nodes = 0;
	FA->dee_clash = 0.5;
	FA->intrafraction = 1.0;
	FA->cluster_rmsd = 2.0f;

	FA->force_interaction=0;
	FA->interaction_factor=5.0;
	FA->atm_cnt=0;
	FA->atm_cnt_real=0;
	FA->res_cnt=0;
	FA->nors=0;
	//FA->natoms_rmsd=0;

    FA->nrg_suite=0;
    FA->nrg_suite_timeout=60;
	FA->translational=0;
	FA->refstructure=0;

	FA->is_protein=1;
	//FA->is_nucleicacid=0;

	FA->delta_angstron=0.25;
	FA->delta_angle=5.0;
	FA->delta_dihedral=5.0;
	FA->delta_flexible=10.0;
	FA->delta_index=1.0;
	FA->max_results=10;

	GB->outgen=0;
	FA->num_grd=0;
	FA->exclude_het=1;
	FA->remove_water=0;

	FA->recalci=0;
	FA->skipped=0;
  
	FA->spacer_length=0.375;
	FA->opt_grid=0;

	FA->pbloops=1;
	FA->bloops=2;

	FA->rotobs=0;

        FA->output_scored_only=0;
	FA->permeability=1.0;
	FA->intramolecular=1;
	FA->solventterm=-2.0;
	FA->by_solventtype=0;
	FA->metaltype=9;        // neutral (AMINO12.def)
    
	FA->useflexdee=0;
	FA->num_constraints=0;

	FA->npar=0;

	FA->mov[0] = NULL;
	FA->mov[1] = NULL;
    
    strcpy(FA->vcontacts_self_consistency,"MEAN");
    
	// Linux path
	pch=strrchr(argv[0],'\\');
	if(pch==NULL) {
		// Windows path
		pch=strrchr(argv[0],'/');
	}

#ifndef _WIN32
	if(pch!=NULL){
		for(i=0;i<(int)(pch-argv[0]);i++){
			FA->base_path[i]=argv[0][i];
			FA->base_path[i+1]='\0';
		}
	}else{
		strcpy(FA->base_path,".");
	}  

#else
	strcpy(FA->base_path,".");
#endif //_WIN32

	printf("base path is '%s'\n", FA->base_path);
  
	strcpy(dockinp,argv[1]);
	strcpy(gainp,argv[2]);
	strcpy(end_strfile,argv[3]);
	strcpy(FA->rrgfile,end_strfile);
	//printf("END FILE:<%s>\n",end_strfile);
	//PAUSE;

	/*
	  if(IS_BIG_ENDIAN())
	  printf("platform is big-endian\n");
	  else
	  printf("platform is little-endian\n");    
	*/

	// use boinc_init();
	Initialize();
  
	wif083(FA);
  
	///////////////////////////////////////////////////////////////////////////////
	// memory allocations for param structures
  
	//printf("memory allocation for opt_par\n");

	FA->map_par = (optmap*)malloc(FA->MIN_PAR*sizeof(optmap));
	FA->opt_par = (double*)malloc(FA->MIN_PAR*sizeof(double));
	FA->del_opt_par = (double*)malloc(FA->MIN_PAR*sizeof(double));
	FA->min_opt_par = (double*)malloc(FA->MIN_PAR*sizeof(double));
	FA->max_opt_par = (double*)malloc(FA->MIN_PAR*sizeof(double));
	FA->map_opt_par = (int*)malloc(FA->MIN_PAR*sizeof(int));

	if(!FA->map_par || !FA->opt_par ||
	   !FA->del_opt_par || !FA->min_opt_par || 
	   !FA->max_opt_par || !FA->map_opt_par){
		fprintf(stderr,"ERROR: memory allocation error for opt_par\n");
		Terminate(2);
	}

	memset(FA->map_par,0,FA->MIN_PAR*sizeof(optmap));
	memset(FA->opt_par,0,FA->MIN_PAR*sizeof(double));
	memset(FA->del_opt_par,0,FA->MIN_PAR*sizeof(double));
	memset(FA->min_opt_par,0,FA->MIN_PAR*sizeof(double));
	memset(FA->max_opt_par,0,FA->MIN_PAR*sizeof(double));

	/////////////////////////////////////////////////////////////////////////////////
  
	printf("Reading input (%s)...\n",dockinp);
	read_input(FA,&atoms,&residue,&rotamer,&cleftgrid,dockinp);
  
	if (strcmp(FA->complf,"VCT")==0){
		// Vcontacts memory allocations...
		// ca_rec can be reallocated
		VC->Calc = (atomsas*)malloc(FA->atm_cnt_real*sizeof(atomsas));
		VC->Calclist = (int*)malloc(FA->atm_cnt_real*sizeof(int));
		VC->ca_index = (int*)malloc(FA->atm_cnt_real*sizeof(int));
		VC->seed = (int*)malloc(3*FA->atm_cnt_real*sizeof(int));
		VC->contlist = (contactlist*)malloc(10000*sizeof(contactlist));
    
		if((!VC->Calc) || (!VC->ca_index) || 
		   (!VC->seed) || (!VC->contlist) || (!VC->Calclist)) {
			fprintf(stderr, "ERROR: memory allocation error for (Calc or Calclist or ca_index or seed or contlist)\n");
			Terminate(2);
		}
	}  
  
	//printf("Create rebuild list...\n");
	create_rebuild_list(FA,atoms,residue);
  
	//printf("atm_cnt=%d\tres_cnt=%d\n",FA->atm_cnt,FA->res_cnt);
	//printf("npar=%d\n",FA->npar);
	//cf=ic2cf(FA,VC,atoms,residue,cleftgrid,FA->npar,FA->opt_par);
	//for(i=0;i<FA->npar;i++){printf("[%8.3f]",FA->opt_par[i]);}
	//printf("=%8.5f\n",cf);

	//-----------------------------------------------------------------------------------
	strcpy(tmp_end_strfile,end_strfile);
	strcat(tmp_end_strfile,"_INI.pdb");
	strcpy(remark,"REMARK initial structure\n");

	VC->first = 1;
	// Should execute cf-vcfunction instead to avoid rotamer change for INI conf.
	cf=ic2cf(FA,VC,atoms,residue,cleftgrid,FA->npar,FA->opt_par);
	VC->first = 0;

	for(i=0;i<FA->npar;i++){printf("[%8.3f]",FA->opt_par[i]);}
	printf("=%8.5f\n", get_cf_evalue(&cf));
	//getchar();
  
	sprintf(tmpremark,"REMARK CF=%8.5f\n", get_cf_evalue(&cf));
	strcat(remark,tmpremark);
	sprintf(tmpremark,"REMARK CF.app=%8.5f\n", get_apparent_cf_evalue(&cf));
	strcat(remark,tmpremark);

	for(i=0;i<FA->num_optres;i++){
    
		res_ptr = &residue[FA->optres[i].rnum];
		cf_ptr = &FA->optres[i].cf;

		sprintf(tmpremark,"REMARK optimizable residue %s %c %d\n",
			res_ptr->name,res_ptr->chn,res_ptr->number);
		strcat(remark,tmpremark);
    
		sprintf(tmpremark,"REMARK CF.com=%8.5f\n",cf_ptr->com);
		strcat(remark,tmpremark);
		sprintf(tmpremark,"REMARK CF.sas=%8.5f\n",cf_ptr->sas);
		strcat(remark,tmpremark);
		sprintf(tmpremark,"REMARK CF.wal=%8.5f\n",cf_ptr->wal);
		strcat(remark,tmpremark);
		sprintf(tmpremark,"REMARK CF.con=%8.5f\n",cf_ptr->con);
		strcat(remark,tmpremark);

	}

	for(i=0;i<FA->npar;i++){
		sprintf(tmpremark,"REMARK [%8.3f]\n",FA->opt_par[i]);
		strcat(remark,tmpremark);
	}
	sprintf(tmpremark,"REMARK inputs: %s & %s\n",dockinp,gainp);
	strcat(remark,tmpremark);
  
	write_pdb(FA,atoms,residue,tmp_end_strfile,remark);

	//printf("wrote initial PDB structure on %s\n",tmp_end_strfile);
	//-----------------------------------------------------------------------------------


	/* PRINTS ALL ACCEPTED ROTAMER LIST
	   for (i=0;i<FA->nflxsc;i++){
	   resnum=FA->flex_res[i].rnum;
	   for (j=1;j<=FA->res_cnt;j++){
	   if (residue[j].number==resnum){
	   printf("ROTAMERS RESIDUE %s%d%c\n-----------------\n",
	   residue[j].name,residue[j].number,residue[j].chn);
	   for (k=0;k<residue[j].trot+1;k++){
	   firstatm=residue[j].fatm[k];
	   lastatm=residue[j].latm[k];
	   printf("Rotamer[%d]\tFATM=%d\tLATM=%d\n",residue[j].rotid[k],firstatm,lastatm);
	   printf("COOR=");
	   for (l=0;l<3;l++){
	   printf("[%1.3f] ",atoms[lastatm].coor[l]);
	   }
	   printf("\n");
	   }
	   }
	   }
	   PAUSE;
	   }
	*/
  
	if(strcmp(FA->metopt,"GA") == 0){

		////////////////////////////////
		////// Genetic Algorithm ///////
		////////////////////////////////

		// calculate time 
		sta_timer=time(NULL);
		sta=localtime(&sta_timer);
		sta_val[0]=sta->tm_sec;
		sta_val[1]=sta->tm_min;
		sta_val[2]=sta->tm_hour;

		int n_chrom_snapshot=GA(FA,GB,VC,&chrom,&chrom_snapshot,&gene_lim,atoms,residue,&cleftgrid,gainp,&memchrom,ic2cf);
    
		if(n_chrom_snapshot > 0){

			end_timer=time(NULL);
			end=localtime(&end_timer);
			end_val[0]=end->tm_sec;
			end_val[1]=end->tm_min;
			end_val[2]=end->tm_hour;
      
			printf("GA:Start time =%0d:%0d:%0d\n",sta_val[2],sta_val[1],sta_val[0]);
			printf("GA:End time   =%0d:%0d:%0d\n",end_val[2],end_val[1],end_val[0]);
      
			ct=0;
			if (sta_val[0]>end_val[0]){
				end_val[1]--;
				end_val[0]+=60;
			}
			if (sta_val[1]>end_val[1]){
				end_val[2]--;
				end_val[1]+=60;
			}
			ct+=((end_val[0]-sta_val[0])+(end_val[1]-sta_val[1])*60);
      
			printf("GA Computational time %ld sec (%4.2f min)\n",ct,(double)ct/60.0);
      
			printf("atoms recalculated=%d\n",FA->recalci);
			printf("individuals skipped=%d\n",FA->skipped);

			////////////////////////////////
			//////       END         ///////
			////////////////////////////////
      
			/******************************************************************/
      
			printf("clustering all individuals in GA...\n");
            fflush(stdout);
            
			cluster(FA,GB,VC,chrom_snapshot,gene_lim,atoms,residue,cleftgrid,n_chrom_snapshot,end_strfile,tmp_end_strfile,dockinp,gainp);
		}
	}
    

	//////////////////////////////////////////
	// free up memory allocated using malloc//
	//////////////////////////////////////////

	// Genes properties
	if(gene_lim != NULL) free(gene_lim);

	// Chromosomes
	if(chrom != NULL){
		for(i=0;i<memchrom;++i){
			if(chrom[i].genes != NULL) free(chrom[i].genes);
		}
		free(chrom);
	}

	if(chrom_snapshot != NULL){
		for(i=0;i<(GB->num_chrom*GB->max_generations);++i){
			if(chrom_snapshot[i].genes != NULL) free(chrom_snapshot[i].genes);
		}
		free(chrom_snapshot);
	}

	// Vcontacts
	if(VC->Calc != NULL) {
		free(VC->Calc);
		free(VC->Calclist);
		free(VC->ca_index);
		free(VC->seed);
		free(VC->contlist);
	}

	// Cleft Grid
	if(cleftgrid != NULL) free(cleftgrid);

	// Atoms
	if(atoms != NULL) {
	  
		for(i=0;i<FA->MIN_NUM_ATOM;i++){

			if(atoms[i].cons != NULL) { free(atoms[i].cons); }
			if(atoms[i].coor_ref != NULL) { free(atoms[i].coor_ref); }

			if(atoms[i].eigen != NULL){
				for(j=0;j<FA->normal_modes;j++)
					if(atoms[i].eigen[j] != NULL)
						free(atoms[i].eigen[j]);
				free(atoms[i].eigen);
			}
		}

		free(atoms);
		
	}
  
	free(FA->num_atm);
  
	for(i=1; i<=FA->ntypes; i++)
		free(FA->energy[i]);
	free(FA->energy);

	// Constraints
	if(FA->constraints != NULL) free(FA->constraints);

	// Residues
	if(residue != NULL) {
		for(i=1;i<=FA->res_cnt;i++){
			//printf("Residue[%d]\n",i);

			if(residue[i].bonded != NULL){
				natm = residue[i].latm[0]-residue[i].fatm[0]+1;
				for(j=0;j<natm;j++){ free(residue[i].bonded[j]); }
				free(residue[i].bonded);
			}
		  
			if(residue[i].gpa != NULL) free(residue[i].gpa);
			if(residue[i].fatm != NULL) free(residue[i].fatm);
			if(residue[i].latm != NULL) free(residue[i].latm);
			if(residue[i].bond != NULL) free(residue[i].bond);
		}

		free(residue);
	}

	// Mov (buildlist)
	for(i=0;i<2;i++){ if(FA->mov[i] != NULL) free(FA->mov[i]); }

	// Optimizable residues
	if(FA->optres != NULL) free(FA->optres);

	// Rotamers
	if(rotamer != NULL) free(rotamer);
  
	// Flexible Residues
	if (FA->flex_res != NULL) {
		for(i=0;i<FA->MIN_FLEX_RESIDUE;i++){
			if(FA->flex_res[i].close != NULL) {
				free(FA->flex_res[i].close);
			}
		}
		free(FA->flex_res);
	}

	// eigen vectors
	if(FA->eigenvector != NULL){
		for(i=0;i<3*FA->MIN_NUM_ATOM;i++)
			if(FA->eigenvector[i] != NULL) 
				free(FA->eigenvector[i]);
		
		free(FA->eigenvector);
	}

	// normal grid
	if(FA->normal_grid != NULL){
		for(i=0;i<FA->MIN_NORMAL_GRID_POINTS;i++)
			if(FA->normal_grid[i] != NULL)
				free(FA->normal_grid[i]);
		
		free(FA->normal_grid);
	}
  
	// Param
	if(FA->map_par != NULL) free(FA->map_par);
	if(FA->opt_par != NULL) free(FA->opt_par);
	if(FA->del_opt_par != NULL) free(FA->del_opt_par);
	if(FA->min_opt_par != NULL) free(FA->min_opt_par);
	if(FA->max_opt_par != NULL) free(FA->max_opt_par);

	/*
	// RMSD
	for(i=0;i<FA->num_het;i++){
	if(FA->res_rmsd[i].fatm != NULL) free(FA->res_rmsd[i].fatm);
	if(FA->res_rmsd[i].latm != NULL) free(FA->res_rmsd[i].fatm);
	}
	//if(FA->atoms_rmsd != NULL) free(FA->atoms_rmsd);
	*/

  
	// FlexDEE Nodes
	if ( FA->psFlexDEENode != NULL ) {
		FA->psFlexDEENode = FA->psFlexDEENode->last;

		while( FA->psFlexDEENode->prev != NULL ) {
			free(FA->psFlexDEENode->rotlist);
			FA->psFlexDEENode = FA->psFlexDEENode->prev;
			free(FA->psFlexDEENode->next);
		}
    
		free(FA->psFlexDEENode->rotlist);
		free(FA->psFlexDEENode);

	}

	if(VC != NULL){
		free(VC->ptorder);
		free(VC->centerpt);
		free(VC->poly);
		free(VC->cont);
		free(VC->vedge);
	  
		free(VC);
	}

	if(GB != NULL) { free(GB); }

	if(FA != NULL) { 
		free(FA->contacts);
		free(FA); 
	}


	//////////////////////////////////////////
	/////////////////  END   /////////////////
	//////////////////////////////////////////

    printf("Done.\n");
    
	Terminate(0);

	return (0);
}
