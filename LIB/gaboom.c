#include "gaboom.h"
#include "Vcontacts.h"
#include "boinc.h"

#ifdef _WIN32
# include <windows.h>
#else
# include <unistd.h>
#endif

using namespace std;

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int GA(FA_Global* FA, GB_Global* GB,VC_Global* VC,chromosome** chrom,chromosome** chrom_snapshot,genlim** gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid,char gainpfile[], int* memchrom, double (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*)){
  
	int i;
	int print=0;
	char tmp_rrgfile[MAX_PATH__];
	char outfile[MAX_PATH__];
	int rrg_flag;
	int rrg_skip=100;
	int n_chrom_snapshot=0;
	char gridfile[MAX_PATH__];

	int geninterval=50;
	int popszpartition=100;
  
	int  state=0;
	char PAUSEFILE[MAX_PATH__];
	char ABORTFILE[MAX_PATH__];
	char STOPFILE[MAX_PATH__];

	const int INTERVAL = 1; // sleep interval between checking file state

	RNGType rng(static_cast<unsigned int>(time(0)));
	boost::uniform_int<> one_to_max_int32( 0, MAX_RANDOM_VALUE );
	boost::variate_generator< RNGType, boost::uniform_int<> >
		dice(rng, one_to_max_int32);

	*memchrom=0; //num chrom allocated in memory

	// for generation random doubles from [0,1[ (mutation crossover operators)
	srand((unsigned)time(0));

	strcpy(PAUSEFILE,FA->state_path);
	strcat(PAUSEFILE,"/.pause");
  
	strcpy(ABORTFILE,FA->state_path);
	strcat(ABORTFILE,"/.abort");

	strcpy(STOPFILE,FA->state_path);
	strcat(STOPFILE,"/.stop");

	GB->num_genes=FA->npar;
	printf("num_genes=%d\n",GB->num_genes);
  
	GB->rrg_skip=0;
	GB->adaptive_ga=0;
	GB->num_print=10;
	GB->print_int=1;
    
    GB->ssnum = 1000;
    GB->pbfrac = 1.0;
    GB->duplicates = 0;
    
	printf("file in GA is <%s>\n",gainpfile);
  
	read_gainputs(FA,GB,gene_lim,&geninterval,&popszpartition,gainpfile);
	printf("read ga inputs\n");
  
	if(GB->print_int < 0){ GB->print_int = 1; }
  
	if(GB->rrg_skip > 0){ rrg_skip = GB->rrg_skip; }
  
	if(GB->num_print > GB->num_chrom){ GB->num_print = GB->num_chrom; }

	if(popszpartition > GB->num_chrom){ popszpartition = GB->num_chrom; }
  
	if(FA->opt_grid){
		printf("will partition grid every %d generations considering %d individuals\n",
		       geninterval, popszpartition);
	}

    double n_poss = set_bins((*gene_lim),GB->num_genes);
    
    printf("%.1lf\n", n_poss);
    if(n_poss < GB->num_chrom && !GB->duplicates){        
        fprintf(stderr,"Too many chromosomes for the number of possibilites (no duplicates allowed).\n");
        fprintf(stderr,"Duplicates are then allowed.\n");
        GB->duplicates = 1;
    }
    
	(*memchrom) = GB->num_chrom;
	if(strcmp(GB->rep_model,"STEADY")==0){
		(*memchrom) += GB->ssnum;
	}else if(strcmp(GB->rep_model,"BOOM")==0){
		(*memchrom) += (int)(GB->pbfrac*(double)GB->num_chrom);
	}

	//printf("memchrom=%d\n",(*memchrom));
	//printf("num_genes=%d\n",GB->num_genes);

	// *** chrom
	(*chrom) = (chromosome*)malloc((*memchrom)*sizeof(chromosome));
	if(!(*chrom)){
		fprintf(stderr,"ERROR: memory allocation error for chrom.\n");
		Terminate(2);
	}
	
	for(i=0;i<(*memchrom);++i){
		(*chrom)[i].genes = (gene*)malloc(GB->num_genes*sizeof(gene));

		if(!(*chrom)[i].genes){
			fprintf(stderr,"ERROR: memory allocation error for chrom[%d].genes.\n",i);	
			Terminate(2);
		}

		(*chrom)[i].evalue = 0.0;
		(*chrom)[i].fitnes = 0.0;
		(*chrom)[i].status = ' ';
	}

	// *** chrom_snapshot
	(*chrom_snapshot) = (chromosome*)malloc((GB->num_chrom*GB->max_generations)*sizeof(chromosome));
	if(!(*chrom_snapshot)){
		fprintf(stderr,"ERROR: memory allocation error for chrom_snapshot.\n");
		Terminate(2);
	}

	for(i=0;i<(GB->num_chrom*GB->max_generations);++i){
		(*chrom_snapshot)[i].genes = (gene*)malloc(GB->num_genes*sizeof(gene));

		if(!(*chrom_snapshot)[i].genes){
			fprintf(stderr,"ERROR: memory allocation error for chrom_snapshot[%d].genes.\n",i);	
			Terminate(2);
		}

		(*chrom_snapshot)[i].evalue = 0.0;
		(*chrom_snapshot)[i].fitnes = 0.0;
		(*chrom_snapshot)[i].status = ' ';
		//printf("chrom_snapshot[%d] allocated at address %p!\n", i, &(*chrom_snapshot)[i]);
	}
    
	printf("alpha %lf peaks %lf scale %lf\n",GB->alpha,GB->peaks,GB->scale);
	GB->sig_share=0.0;
  
	for(i=0;i<GB->num_genes;i++){
		//printf("max=%6.3f min=%6.3f del=%6.3f\n",(*gene_lim)[i].max,(*gene_lim)[i].min,(*gene_lim)[i].del);
		//PAUSE;
		GB->sig_share += ((*gene_lim)[i].max-(*gene_lim)[i].min)*((*gene_lim)[i].max-(*gene_lim)[i].min);
	}
	GB->sig_share = sqrt(GB->sig_share/(double)GB->num_genes)/(2.0*pow(GB->peaks,(1.0/(double)GB->num_genes)));
	GB->sig_share /= GB->scale;
	printf("SIGMA_SHARE=%f\n",GB->sig_share);
  
  
	// for(i=0;i<GB->num_genes;i++) {
	//printf("GENE(%d)=[%8.3f,%8.3f,%8.3f,%d]\n",
	//	   i,(*gene_lim)[i].min,(*gene_lim)[i].max,(*gene_lim)[i].del);
	//PAUSE;
  
	populate_chromosomes(FA,GB,VC,(*chrom),(*gene_lim),atoms,residue,(*cleftgrid),GB->pop_init_method,target,GB->pop_init_file,0,print,dice);
	//}
	
	//print_pop((*chrom),(*gene_lim),GB->num_chrom,GB->num_genes);
  
	/*
	  for(i=0;i<GB->num_genes;i++){
	  printf("%d %f %f %f\n",i,GB->min_opt_par[i],GB->max_opt_par[i],GB->del_opt_par[i]);
	  PAUSE;
	  }
	*/
    
        int save_num_chrom = (int)(GB->num_chrom*SAVE_CHROM_FRACTION);

	for(i=0;i<GB->max_generations;i++){

		///////////////////////////////////////////////////

		state=check_state(PAUSEFILE,ABORTFILE,STOPFILE,INTERVAL);
    
		if(state == -1){ 
			return(state); 

		}else if(state == 1){ 
			
			quicksort_evalue((*chrom_snapshot),0,n_chrom_snapshot-1);
			return(n_chrom_snapshot); 
		}

		////////////////////////////////

		// BOINC CLIENT GUI PROGRESS BAR
#ifdef ENABLE_BOINC
		boinc_fraction_done((double)(i+1)/(double)GB->max_generations);
#endif

		////////////////////////////////
		
		//printf("chrom_snapshot[%d] at address %p\n", i*GB->num_chrom, chrom_snapshot[i*GB->num_chrom]);
		if (FA->opt_grid                    &&     // if a OPTGRD line was specified
		    ((i+1) % geninterval) == 0      &&     // is factor of
		    (i+1) != GB->max_generations) {        // discard the last generation
      
			//need to sort in decreasing order of energy
			quicksort_evalue((*chrom),0,GB->num_chrom-1);
      
			//printf("Partionning grid...(%d)\n",FA->popszpartition);
			partition_grid(FA,(*chrom),(*gene_lim),atoms,residue,cleftgrid,popszpartition,1);
      
			if(FA->output_range){
				sprintf(gridfile,"grid.%d.prt.pdb",i+1);
				write_grid(FA,(*cleftgrid),gridfile);
			}
      
			slice_grid(FA,(*gene_lim),atoms,residue,cleftgrid);
      
			if(FA->output_range){
				sprintf(gridfile,"grid.%d.slc.pdb",i+1);
				write_grid(FA,(*cleftgrid),gridfile);
			}
      
			//repopulate unselected individuals
			populate_chromosomes(FA,GB,VC,(*chrom),(*gene_lim),atoms,residue,(*cleftgrid),GB->pop_init_method,target,GB->pop_init_file,popszpartition,print,dice);
		}
    
		print = ( (i+1) % GB->print_int == 0 ) ? 1 : 0;
		if(print) { printf("Generation: %3d\n",i+1); }

		//print_par(chrom,gene_lim,20,GB->num_genes);
		//PAUSE;

		rrg_flag=0;
		if((i/rrg_skip)*rrg_skip == i) rrg_flag=1;
		if((rrg_flag==1) && (GB->outgen==1)){
			if(FA->refstructure == 1){
				sprintf(tmp_rrgfile,"%s_%d.rrg",FA->rrgfile,i);
				//printf("%s\n",tmp_rrgfile);
				//PAUSE;
				write_rrg(FA,GB,(*chrom),(*gene_lim),atoms,residue,(*cleftgrid),tmp_rrgfile);
			}
			strcpy(outfile,"gaboom_par.res");
			write_par((*chrom),(*gene_lim),i,outfile,GB->num_chrom,GB->num_genes);
		}

		reproduce(FA,GB,VC,(*chrom),(*gene_lim),atoms,residue,(*cleftgrid),GB->rep_model,GB->mut_rate,GB->cross_rate,print,dice,target);
		
		save_snapshot(&(*chrom_snapshot)[i*GB->num_chrom],(*chrom),save_num_chrom,GB->num_genes);
		n_chrom_snapshot += save_num_chrom;
		

		if(strcmp(GB->fitness_model,"PSHARE")==0){
			quicksort_fitnes((*chrom),0,GB->num_chrom-1);
			
			if(print){
				printf("best by fitnes\n");
				print_par((*chrom),(*gene_lim),GB->num_print,GB->num_genes);
			}
		}
	}
 
	quicksort_evalue((*chrom),0,GB->num_chrom-1);

	if (GB->outgen) write_par((*chrom),(*gene_lim),GB->max_generations,outfile,GB->num_chrom,GB->num_genes);

	quicksort_evalue((*chrom_snapshot),0,n_chrom_snapshot-1);

	/*
		printf("Save snapshot == END ==\n");
		print_par((*chrom_snapshot),(*gene_lim),n_chrom_snapshot,GB->num_genes);
	*/

	n_chrom_snapshot = remove_dups((*chrom_snapshot),n_chrom_snapshot,GB->num_genes);
	
	/*	
		printf("Save snapshot == END ==\n");
		print_par((*chrom_snapshot),(*gene_lim),n_chrom_snapshot,GB->num_genes);
	*/
	
	return n_chrom_snapshot;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void save_snapshot(chromosome* chrom_snapshot, const chromosome* chrom, int num_chrom, int num_genes){
	
	for(int i=0; i<num_chrom; i++){
		chrom_snapshot[i].evalue = chrom[i].evalue;
		chrom_snapshot[i].fitnes = chrom[i].fitnes;
		chrom_snapshot[i].status = chrom[i].status;
		
		for(int j=0; j<num_genes; j++){
			chrom_snapshot[i].genes[j].to_ic = chrom[i].genes[j].to_ic;
			chrom_snapshot[i].genes[j].to_int32 = chrom[i].genes[j].to_int32;
		}
		
	}

}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int check_state(char* pausefile, char* abortfile, char* stopfile, int interval){
	FILE* STATE;
  
	STATE = NULL;

	// try and open pause/stop file
	// (works with PyMOL interface)
  
	STATE = fopen(pausefile,"r");
	if(STATE != NULL) {
		do {
			fclose(STATE);
			STATE = fopen(pausefile,"r");

# ifdef _WIN32
			Sleep(1000);
# else
			sleep(1);
# endif

		}while(STATE != NULL);
	}
  
	STATE = fopen(abortfile,"r");
	if(STATE != NULL) {
		fclose(STATE);
		printf("manual aborting\n");
		return -1;
	}

	STATE = fopen(stopfile,"r");
	if(STATE != NULL) {
		fclose(STATE);
		printf("simulation stopped prematurely\n");
		return 1;
	}
  
	return 0;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void fitness_stats(GB_Global* GB, const chromosome* chrom,int pop_size){
	int i;
	int flag;

	//calculate fitness max and and average of the whole pop
	GB->fit_max=0.0;
	GB->fit_avg=0.0;

	flag=1;
	for(i=0;i<pop_size-i;i++){
		if (flag){
			GB->fit_max=chrom[i].fitnes;
			flag=0;
		}

		if(chrom[i].fitnes > GB->fit_max) 
			GB->fit_max=chrom[i].fitnes;
    
		GB->fit_avg+=chrom[i].fitnes;
	}
  
	GB->fit_avg /= (double)pop_size;
  
	return;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void adapt_prob(GB_Global* GB,double fit1, double fit2, double* mutp, double* crossp){
	//printf("crossing fit1[%8.3f] with fit2[%8.3f]\n",fit1,fit2);
  
	//find which crossed individual has higher fitness
	if(fit1 > fit2){
		GB->fit_high=fit1;
		GB->fit_low=fit2;
	}else{
		GB->fit_high=fit2;
		GB->fit_low=fit1;
	}

	//crossp=k1 when high=avg
	//mutp=k2 when high=avg
	//crossp/mutp=0 when high=max

	//calculate new probabilities (pc/pm)
	if (GB->fit_high > GB->fit_avg) *crossp = GB->k1*(GB->fit_max-GB->fit_high)/(GB->fit_max-GB->fit_avg);
	else *crossp = GB->k3;
  
	if (GB->fit_low > GB->fit_avg) *mutp = GB->k2*(GB->fit_max-GB->fit_low)/(GB->fit_max-GB->fit_avg);
	else *mutp = GB->k4;
  
	/*
	  printf("f'=%.1f\tf=%.1f\tfmax=%.1f\tfavg=%.1f\t\tPc=%5.3f\tPm=%5.3f\n",
	  GB->fit_high,GB->fit_low,GB->fit_max,GB->fit_avg,*crossp,*mutp);
	*/

	return;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void reproduce(FA_Global* FA,GB_Global* GB,VC_Global* VC, chromosome* chrom,const genlim* gene_lim,
               atom* atoms,resid* residue,gridpoint* cleftgrid,char* repmodel, 
               double mutprob, double crossprob, int print, 
               boost::variate_generator< RNGType, boost::uniform_int<> > &, 
               double (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*)){

	int i,j,k;
	int nnew,p1,p2;

	gene chrop1_gen[MAX_NUM_GENES];
	gene chrop2_gen[MAX_NUM_GENES];

	int num_genes_wo_sc=0;	

	//before calculating get avg and max fitness of the whole pop.
	if (GB->adaptive_ga) {
		fitness_stats(GB,chrom,GB->num_chrom);
		/*
		  printf("------fitness stats-------\navg=%8.3f\tmax=%8.3f\n",GB->fit_avg,GB->fit_max);
		  printf("k1=%5.3f\tk2=%5.3f\tk3=%5.3f\tk4=%5.3f\n",GB->k1,GB->k2,GB->k3,GB->k4);
		  PAUSE;
		*/
	}

    //***************************ELITISM**********************************
    
	if(strcmp(repmodel,"STEADY")==0){

		i=0;
		while(i<GB->ssnum){

			/************************************/
			/****** SELECTION OF PARENTS ********/
			/************************************/

			p1=roullete_wheel(chrom,GB->num_chrom);
			p2=roullete_wheel(chrom,GB->num_chrom);
			if (GB->adaptive_ga) adapt_prob(GB,chrom[p1].fitnes,chrom[p2].fitnes,&mutprob,&crossprob);

			/************************************/
			/******  CROSSOVER OPERATOR  ********/
			/************************************/

			// create temporary genes
			memcpy(chrop1_gen,chrom[p1].genes,GB->num_genes*sizeof(gene));
			memcpy(chrop2_gen,chrom[p2].genes,GB->num_genes*sizeof(gene));

			if(RandomDouble() < crossprob){
				crossover(chrop1_gen,chrop2_gen,GB->num_genes);
			}

			/************************************/
			/******   MUTATION OPERATOR  ********/
			/************************************/

			num_genes_wo_sc = GB->num_genes-FA->nflxsc_real;

			mutate(chrop1_gen,GB->num_genes-FA->nflxsc_real,mutprob);
			k=0;
			for(j=0;j<FA->nflxsc;j++){
				if(residue[FA->flex_res[j].inum].trot != 0){
					if(RandomDouble() < FA->flex_res[j].prob){
						mutate(&chrop1_gen[num_genes_wo_sc+k],1,mutprob);
					}
					k++;
				}
			}

			mutate(chrop2_gen,GB->num_genes-FA->nflxsc_real,mutprob);
			k=0;
			for(j=0;j<FA->nflxsc;j++){
				if(residue[FA->flex_res[j].inum].trot != 0){
					if(RandomDouble() < FA->flex_res[j].prob){
						mutate(&chrop1_gen[num_genes_wo_sc+k],1,mutprob);
					}
					k++;
				}
			}

			for(j=0; j<GB->num_genes; j++){
				chrop1_gen[j].to_ic = genetoic(&gene_lim[j],chrop1_gen[j].to_int32);
				chrop2_gen[j].to_ic = genetoic(&gene_lim[j],chrop2_gen[j].to_int32);
			}

			/************************************/
			/******   CHECK DUPLICATION  ********/
			/************************************/

			if(GB->duplicates || cmp_chrom2pop(chrom,chrop1_gen,GB->num_genes,0,GB->num_chrom+i)==0){

				if(!FA->useflexdee || cmp_chrom2rotlist(FA->psFlexDEENode,chrom,gene_lim,num_genes_wo_sc,FA->nflxsc_real,GB->num_chrom,FA->FlexDEE_Nodes)==0){
	  
					memcpy(chrom[GB->num_chrom+i].genes,chrop1_gen,GB->num_genes*sizeof(gene));
					chrom[GB->num_chrom+i].evalue=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,chrom[GB->num_chrom+i].genes,target);
					chrom[GB->num_chrom+i].status='n';
					i++;
	  
				}
			}
      
			if(i==GB->ssnum) break;
      
			if(GB->duplicates || cmp_chrom2pop(chrom,chrop2_gen,GB->num_genes,0,GB->num_chrom+i)==0){
	
				if(!FA->useflexdee || cmp_chrom2rotlist(FA->psFlexDEENode,chrom,gene_lim,num_genes_wo_sc,FA->nflxsc_real,GB->num_chrom,FA->FlexDEE_Nodes)==0){
	  
					memcpy(chrom[GB->num_chrom+i].genes,chrop2_gen,GB->num_genes*sizeof(gene));
					chrom[GB->num_chrom+i].evalue=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,chrom[GB->num_chrom+i+1].genes,target);
					chrom[GB->num_chrom+i].status='n';
					i++;
	  
				}
			}
		}
		for(i=0;i<GB->ssnum;i++) chrom[GB->num_chrom-1-i]=chrom[GB->num_chrom+i];
		calculate_fitness(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,GB->fitness_model,GB->num_chrom,print,target);
	}

	//****************** POPULATION BOOM ****************

	if(strcmp(repmodel,"BOOM")==0){
    
		/* create n new individuals from old pop according to fitness
		   without duplicates, sort them together with the whole pop
		   of N+n individuals and select the N best fit */
		nnew=(int)(GB->pbfrac*(double)GB->num_chrom);
		i=0;
		while(i<nnew){

			/************************************/
			/****** SELECTION OF PARENTS ********/
			/************************************/

			p1=roullete_wheel(chrom,GB->num_chrom);
			p2=roullete_wheel(chrom,GB->num_chrom);
			if (GB->adaptive_ga) { adapt_prob(GB,chrom[p1].fitnes,chrom[p2].fitnes,&mutprob,&crossprob); }

			/************************************/
			/******  CROSSOVER OPERATOR  ********/
			/************************************/

			memcpy(chrop1_gen,chrom[p1].genes,GB->num_genes*sizeof(gene));
			memcpy(chrop2_gen,chrom[p2].genes,GB->num_genes*sizeof(gene));

			if(RandomDouble() < crossprob) {
				crossover(chrop1_gen,chrop2_gen,GB->num_genes);
			}

			/************************************/
			/******   MUTATION OPERATOR  ********/
			/************************************/

			num_genes_wo_sc = GB->num_genes-FA->nflxsc_real;

			mutate(chrop1_gen,GB->num_genes-FA->nflxsc_real,mutprob);
			//printf("mutating %p to %d genes\n",chrop1_gen,GB->num_genes-FA->nflxsc_real);

			k=0;
			for(j=0;j<FA->nflxsc;j++){
				if(residue[FA->flex_res[j].inum].trot > 0){
					if(RandomDouble() < FA->flex_res[j].prob){
						//printf("mutating residue %s %c %d with address %p\n", FA->flex_res[j].name, FA->flex_res[j].chn, FA->flex_res[j].num,&chrop1_gen[num_genes_wo_sc+k]);
						mutate(&chrop1_gen[num_genes_wo_sc+k],1,mutprob);
					}
					k++;
				}
			}

			mutate(chrop2_gen,GB->num_genes-FA->nflxsc_real,mutprob);
			k=0;
			for(j=0;j<FA->nflxsc;j++){
				if(residue[FA->flex_res[j].inum].trot > 0){
					if(RandomDouble() < FA->flex_res[j].prob)
						mutate(&chrop2_gen[num_genes_wo_sc+k],1,mutprob);
					k++;
				}
			}

			for(j=0; j<GB->num_genes; j++){
				chrop1_gen[j].to_ic = genetoic(&gene_lim[j],chrop1_gen[j].to_int32);
				chrop2_gen[j].to_ic = genetoic(&gene_lim[j],chrop2_gen[j].to_int32);
			}

			/************************************/
			/******   CHECK DUPLICATION  ********/
			/************************************/
      
			if(GB->duplicates || cmp_chrom2pop(chrom,chrop1_gen,GB->num_genes,0,GB->num_chrom+i)==0){
	
				if(!FA->useflexdee || 
				   cmp_chrom2rotlist(FA->psFlexDEENode,chrom,gene_lim,num_genes_wo_sc,FA->nflxsc_real,GB->num_chrom,FA->FlexDEE_Nodes)==0){

					memcpy(chrom[GB->num_chrom+i].genes,chrop1_gen,GB->num_genes*sizeof(gene));
	  
					chrom[GB->num_chrom+i].evalue=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,chrom[GB->num_chrom+i].genes,target);
					chrom[GB->num_chrom+i].status='n';
					i++;
	  
				}
			}

			if(i==nnew) break;

			if(GB->duplicates || cmp_chrom2pop(chrom,chrop2_gen,GB->num_genes,0,GB->num_chrom+i)==0){

				if(!FA->useflexdee || cmp_chrom2rotlist(FA->psFlexDEENode,chrom,gene_lim,num_genes_wo_sc,FA->nflxsc_real,GB->num_chrom,FA->FlexDEE_Nodes)==0){
	  
					memcpy(chrom[GB->num_chrom+i].genes,chrop2_gen,GB->num_genes*sizeof(gene));
	  
					chrom[GB->num_chrom+i].evalue=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,chrom[GB->num_chrom+i].genes,target);
					chrom[GB->num_chrom+i].status='n';
					i++;

				}
			}
		}
    
		calculate_fitness(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,GB->fitness_model,GB->num_chrom+nnew,print,target);
    	}

	return;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int roullete_wheel(const chromosome* chrom,int n){
	double r;
	double tot=0.0;
	int i;

	for(i=0;i<n;i++){tot += chrom[i].fitnes;}
	//printf("tot=%f\n",tot);
	//PAUSE;

	r=RandomDouble()*tot;
  
	i=0;
	tot=0.0;
	while(tot <= r){
		tot += chrom[i].fitnes;
		i++;
	}
	//printf("r=%f tot=%f i=%d\n",r,tot,i);
	i--;

	//PAUSE;
  
	return i;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void calculate_fitness(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim,
                       atom* atoms,resid* residue,gridpoint* cleftgrid,char method[], int pop_size, int print,
                       double (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*)){

	int i,j;
	//float tot=0.0;
	double share,rmsp;

	for(i=0;i<pop_size;i++){
		if(chrom[i].status != 'n'){
			chrom[i].evalue=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,chrom[i].genes,target);
			chrom[i].status='n';
		}
	}

	quicksort_evalue(chrom,0,pop_size-1);

	//print_par(chrom,gene_lim,5,GB->num_genes);
	//PAUSE;
	//chrom_hpsort(pop_size,0,chrom);
  
	if(strcmp(method,"LINEAR")==0){
		/* the fitness value is a number between 0 and num_chrom.
		   each chromosome is assigned an integer value that 
		   corresponds to its position in index_map.
		*/
		for(i=0;i<GB->num_chrom;i++){
			chrom[i].fitnes=(double)(GB->num_chrom-i);
		}
	}
  
	if(strcmp(method,"PSHARE")==0){
		/* the fitness value is a number between 0 and num_chrom.
		   each chromosome is assigned an integer value that 
		   corresponds to its position in index_map. Moreover,
		   each chromosome's fitness is lowered by sharing.
		*/
    
		for(i=0;i<GB->num_chrom;i++){
      
			share=0.0;
			for(j=0;j<GB->num_chrom;j++){
				
				rmsp=calc_rmsp(GB->num_genes,chrom[i].genes,chrom[j].genes);
				//printf("i=%d j=%d rmsp=%f\n",i,j,rmsp);
				if(rmsp <= GB->sig_share){
					share += (1.0 - pow((rmsp/GB->sig_share),GB->alpha));
				}
				//share=1.0;
				chrom[i].fitnes = (GB->num_chrom-i)/share;
				//printf("i=%d lf=%d share=%f fit=%f\n",i,(GB->num_chrom-i),
				//     share,chrom[i].fitnes);
				//PAUSE;
			}
		}
	}

	if(print){
		printf("best by energy\n");
		print_par(chrom,gene_lim,GB->num_print,GB->num_genes);
	}

	return;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
double eval_chromosome(FA_Global* FA,GB_Global* GB,VC_Global* VC,const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,gene* john, double (*function)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*)){
	
	double icv[MAX_NUM_GENES] = {0};
	
	for(int i=0;i<GB->num_genes;i++){
		if(john[i].to_ic > gene_lim[i].max) {
			fprintf(stderr, "Exceptional out of bounds error at: max: %.5lf when ic: %.5lf\n", gene_lim[i].max, john[i].to_ic);
			john[i].to_ic = gene_lim[i].max;
		}else if(john[i].to_ic < gene_lim[i].min) {
			fprintf(stderr, "Exceptional out of bounds error at: min: %.5lf when ic: %.5lf\n", gene_lim[i].max, john[i].to_ic);
			john[i].to_ic = gene_lim[i].min;
		}

		icv[i] = john[i].to_ic;
	}
  
	return (*function)(FA,VC,atoms,residue,cleftgrid,GB->num_genes,icv);
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void populate_chromosomes(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim, 
                          atom* atoms,resid* residue,gridpoint* cleftgrid,char method[], 
                          double (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*), 
                          char file[], int popoffset, int print, 
                          boost::variate_generator< RNGType, boost::uniform_int<> > &){
  
	int i,j,l;
	
	
	RNGType rng;
	
	boost::uniform_int<> one_to_max_int32( 1, MAX_RANDOM_VALUE );
	boost::variate_generator< RNGType, boost::uniform_int<> >
		dice(rng, one_to_max_int32);
	

	/*
	char *buffer;            
	int size_of_buffer;

	buffer = (char*)malloc((GB->num_genes*13+45)*sizeof(char)); 
	if(!buffer){
		fprintf(stderr,"ERROR: could not allocate memory for buffer\n");
		Terminate(2);
	}
	size_of_buffer= (GB->num_genes*13+45)*sizeof(char);
	*/

	// initialise genes to zero
	for(i=popoffset;i<GB->num_chrom;i++){
		for(j=0;j<GB->num_genes;j++){
			chrom[i].genes[j].to_int32=0;
			chrom[i].genes[j].to_ic=0.0;
		}
	}

	//------------------------------------------------------------------------------
	// use method to create new genes
	if(strcmp(method,"RANDOM")==0){
		printf("generating random population...\n");
		//printf("num_chrom=%d num_genes=%d\n",GB->num_chrom,GB->num_genes);
    
		i=popoffset;
		while(i<GB->num_chrom){
			for(j=0;j<GB->num_genes;j++){

				// side-chain optimization
				if(FA->map_par[j].typ == 4){
					l=0;
					while(FA->flex_res[l].inum != atoms[FA->map_par[j].atm].ofres){
						l++;
					};
	  
					//printf("probability of atom[%d].ofres[%d]\t flex_res[%d](%s).inum[%d]= %.3f\n", FA->map_par[j].atm, atoms[FA->map_par[j].atm].ofres, l, FA->flex_res[l].name, FA->flex_res[l].inum, FA->flex_res[l].prob);

					if(RandomDouble() < FA->flex_res[l].prob){
						chrom[i].genes[j].to_int32 = dice();
					}else{
						chrom[i].genes[j].to_int32 = 0;
					}
				}else{
					chrom[i].genes[j].to_int32 = dice();
					//printf("chrom[%d].gene[%d]=%u\n",i,j,chrom[i].genes[j]);
				}
				
				chrom[i].genes[j].to_ic = genetoic(&gene_lim[j],chrom[i].genes[j].to_int32);
			}
			
			if(GB->duplicates || cmp_chrom2pop(chrom,chrom[i].genes,GB->num_genes,0,i)==0){i++;}
		}
	}

	//------------------------------------------------------------------------------
	/*
	if(strcmp(method,"IPFILE")==0){
		infile_ptr=NULL;
		if(!OpenFile_B(file,"r",&infile_ptr))
			Terminate(8);
    
		i=0;
		while (fgets(buffer,size_of_buffer,infile_ptr)){
			//printf("--->%s<---",buffer);
			//PAUSE;
			if(buffer[3] == '('){
				for(j=0;j<GB->num_genes;j++){
					strncpy(buf_gene,buffer+j*13+4,13);
					buf_gene[13]='\0';
					sscanf(buf_gene,"%f",&gene_float);	  
					chrom[i].genesrand[j]= (int)(
						((gene_float-gene_lim[j].min)/(gene_lim[j].max-gene_lim[j].min))*
						(pow(2.0,(double)gene_blen[j])-1.0));
				}
				i++;      
			}
		}
    
		while(i<GB->num_chrom){
			for(j=0;j<GB->num_genes;j++){
				chrom[i].genes[j]= (int)(RandomDouble()*pow(2.0,(double)gene_blen[j]));
				//printf("chrom[%d].gene[%d]=%d\n",i,j,chrom[i].genes[j]);
			}
			if(GB->duplicates || cmp_chrom2pop(chrom,chrom[i].genes,GB->num_genes,0,GB->num_chrom+i)==0){i++;}
		}
    
		CloseFile_B(&infile_ptr,"r");
	}
	*/
	//------------------------------------------------------------------------------

	// calculate evalue for each chromosome
	for(i=popoffset;i<GB->num_chrom;i++){
		chrom[i].evalue=eval_chromosome(FA,GB,VC,gene_lim,atoms,residue,cleftgrid,chrom[i].genes,target);
		chrom[i].status='n';
		//PAUSE;
		//printf("evalue(%d)=%6.3f\n",i,chrom[i].evalue);
	}
  
	// sort and calculate fitness
	calculate_fitness(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,GB->fitness_model,GB->num_chrom,print,target);

	//free(buffer);
  
	return;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int cmp_chrom2rotlist(psFlexDEE_Node psFlexDEE_INI_Node, const chromosome* chrom, const genlim* gene_lim,int gene_offset, int num_genes, int tot, int num_nodes){

	int   par[100];  
	//int* genes = NULL;
	sFlexDEE_Node sFlexDEENode;

	memset(&par,0,sizeof(par));
  
	if ( psFlexDEE_INI_Node == NULL ) { return 0; }
  

	sFlexDEENode.rotlist = par;

	for(int i=0;i<tot;i++){
		//genes = &chrom[i].genesic[gene_offset];

		psFlexDEE_INI_Node = psFlexDEE_INI_Node->last;
    
		if ( dee_pivot(&sFlexDEENode,&psFlexDEE_INI_Node,1,num_nodes,(num_nodes+1)/2,num_nodes,num_genes) == 0 ) { return 1; }
    
	}
  
	return 0;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int cmp_chrom2pop(const chromosome* chrom,const gene* genes, int num_genes,int start, int last){
	int i,j,flag;
  
	for(i=start;i<last;i++){
		flag=0;
		for(j=0;j<num_genes;j++){
			//printf("comparing %u to %u\n",c->genes[j],chrom[i].genes[j]);
			flag += abs(genes[j].to_ic - chrom[i].genes[j].to_ic) < 0.1;
		}
    
		//printf("flag=%d\n",flag);
		if(flag == num_genes){return 1;}
	}

	return 0;
} 
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int cmp_chrom2pop_int(const chromosome* chrom,const gene* genes, int num_genes,int start, int last){
	int i,j,flag;
  
	for(i=start;i<last;i++){
		flag=0;
		for(j=0;j<num_genes;j++){
			//printf("comparing %u to %u\n",c->genes[j],chrom[i].genes[j]);
			flag += genes[j].to_int32 == chrom[i].genes[j].to_int32;
		}
    
		//printf("flag=%d\n",flag);
		if(flag == num_genes){return 1;}
	}

	return 0;
} 
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
double set_bins(genlim* gene_lim, int num_genes){	
	
    double n_poss = 0.0;
    
	for(int i=0; i<num_genes; i++){
		double nbin = (gene_lim[i].max - gene_lim[i].min) / gene_lim[i].del;
		if(nbin - (int)nbin > 0.0){ nbin += 1.0; }
		if(gene_lim[i].map){ nbin += 1.0; }
		
		gene_lim[i].bin = 1.0/nbin;
            
        if(n_poss > 0.0){
            n_poss *= nbin;
        }else{
            n_poss = nbin;
        }
	}

	return n_poss;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void set_bins(genlim* gene_lim){	
	
	double nbin = (gene_lim->max - gene_lim->min) / gene_lim->del;
	if(nbin - (int)nbin > 0.0){ nbin += 1.0; }
	if(gene_lim->map){ nbin += 1.0; }
	
	gene_lim->bin = 1.0/nbin;
	
	return;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void read_gainputs(FA_Global* FA,GB_Global* GB,genlim** gene_lim,int* gen_int,int* sz_part,char file[]){

	FILE *infile_ptr;        /* pointer to input file */
	char buffer[81];         /* a line from the INPUT file */
	char field[9];           /* field names on INPUT file */
	int ngenes=0;
	int i;

	//printf("file here is <%s>\n",file);
	infile_ptr=NULL;
	if(!OpenFile_B(file,"r",&infile_ptr))
		Terminate(8);
  
	while (fgets(buffer, sizeof(buffer),infile_ptr)){
		for(i=0;i<=7;i++){
			field[i]=buffer[i];
		}
		field[8]='\0';
    
		if(strcmp(field,"NUMCHROM") == 0){
			sscanf(buffer,"%s %d",field,&GB->num_chrom);
		}else if(strcmp(field,"OPTIGRID")==0){
			sscanf(buffer,"%s %d %d %d",field,&FA->opt_grid,gen_int,sz_part);
		}else if(strcmp(field,"NUMGENER") == 0){
			sscanf(buffer,"%s %d",field,&GB->max_generations);
		}else if(strcmp(field,"ADAPTVGA") == 0){
			sscanf(buffer,"%s %d",field,&GB->adaptive_ga);
		}else if(strcmp(field,"ADAPTKCO") == 0){
			//adaptive response parameters
			//k1-k4 are values ranging from 0.0-1.0 inclusively
			sscanf(buffer,"%s %lf %lf %lf %lf",field,&GB->k1,&GB->k2,&GB->k3,&GB->k4);
		}else if(strcmp(field,"CROSRATE") == 0){
			sscanf(buffer,"%s %lf",field,&GB->cross_rate);
		}else if(strcmp(field,"MUTARATE") == 0){
			sscanf(buffer,"%s %lf",field,&GB->mut_rate);
		}else if(strcmp(field,"INIMPROB") == 0){
			sscanf(buffer,"%s %lf",field,&GB->ini_mut_prob);
		}else if(strcmp(field,"ENDMPROB") == 0){
			sscanf(buffer,"%s %lf",field,&GB->end_mut_prob);
		}else if(strcmp(field,"POPINIMT") == 0){      
			sscanf(buffer,"%s %s %s",field,GB->pop_init_method,GB->pop_init_file);
			//printf("<%s> <%s>\n",GB->pop_init_method,GB->pop_init_file);
			//PAUSE;
		}else if(strcmp(field,"FITMODEL") == 0){
			sscanf(buffer,"%s %s",field,GB->fitness_model);
		}else if(strcmp(field,"REPMODEL") == 0){
			sscanf(buffer,"%s %s",field,GB->rep_model);
		}else if(strcmp(field,"DUPLICAT") == 0){
            GB->duplicates = 1;
		}else if(strcmp(field,"BOOMFRAC") == 0){
			sscanf(buffer,"%s %lf",field,&GB->pbfrac);
		}else if(strcmp(field,"STEADNUM") == 0){
			sscanf(buffer,"%s %d",field,&GB->ssnum);
		}else if(strcmp(field,"SHAREALF") == 0){
			sscanf(buffer,"%s %lf",field,&GB->alpha);
		}else if(strcmp(field,"SHAREPEK") == 0){
			sscanf(buffer,"%s %lf",field,&GB->peaks);
		}else if(strcmp(field,"SHARESCL") == 0){
			sscanf(buffer,"%s %lf",field,&GB->scale);
		}else if(strcmp(field,"OUTGENER") == 0){
			GB->outgen = 1;
		}else if(strcmp(field,"PRINTCHR") == 0){
			sscanf(buffer,"%s %d",field,&GB->num_print);
		}else if(strcmp(field,"PRINTINT") == 0){
			sscanf(buffer,"%s %d",field,&GB->print_int);
		}else if(strcmp(field,"PRINTRRG") == 0){
			sscanf(buffer,"%s %d",field,&GB->rrg_skip);
		}else{
			// ...
		}

	}
	CloseFile_B(&infile_ptr,"r");

  
	(*gene_lim) = (genlim*)malloc(GB->num_genes*sizeof(genlim));
	if(!(*gene_lim)){
		fprintf(stderr,"ERROR: memory allocation error for gene_lim.\n");	
		Terminate(2);
	}

	for(ngenes=0;ngenes<GB->num_genes;ngenes++){
		(*gene_lim)[ngenes].min=FA->min_opt_par[ngenes];
		(*gene_lim)[ngenes].max=FA->max_opt_par[ngenes];
		(*gene_lim)[ngenes].del=FA->del_opt_par[ngenes];
		(*gene_lim)[ngenes].map=FA->map_opt_par[ngenes];

		printf("gene %d: min: %8.2f max: %8.2f delta: %8.2f map: %d\n", ngenes,
		       (*gene_lim)[ngenes].min,
		       (*gene_lim)[ngenes].max,
		       (*gene_lim)[ngenes].del,
		       (*gene_lim)[ngenes].map);
    
	}

	return;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void crossover(gene *john,gene *mary,int num_genes){

	/* john and mary are two chromosomes to be crossover at points a and b
	 */

	/*
	  printf("john before:\n");
	  print_chrom(john,num_genes,0);
	  printf("mary before:\n");
	  print_chrom(mary,num_genes,0);
	*/

	int i,j;
	int optr,temp;
	int gen_a,gen_b,pnt_a,pnt_b,aux_gen,aux_pnt;

	gen_a=(int)(RandomDouble()*(double)num_genes);
	pnt_a=(int)(RandomDouble()*(double)(MAX_GEN_LENGTH));
	gen_b=(int)(RandomDouble()*(double)num_genes);
	pnt_b=(int)(RandomDouble()*(double)(MAX_GEN_LENGTH));
	//printf("gen_a=%d\tpnt_a=%d\tgen_b=%d\tpnt_b=%d\n",gen_a,pnt_a,gen_b,pnt_b);

	if(gen_a > gen_b){
		aux_gen=gen_a;
		aux_pnt=pnt_a;
		gen_a=gen_b;
		pnt_a=pnt_b;
		gen_b=aux_gen;
		pnt_b=aux_pnt;
	}

	if(gen_a == gen_b && pnt_a < pnt_b){
		aux_pnt=pnt_b;
		pnt_b=pnt_a;
		pnt_a=aux_pnt;
	}

	for(j=gen_a;j<=gen_b;j++){
		optr=1;
		aux_pnt = (j==gen_a)?pnt_a:(MAX_GEN_LENGTH);
		for(i=0;i<aux_pnt;i++) optr |= (optr << 1);
		temp = (john[j].to_int32 & ~optr) | (mary[j].to_int32 &  optr);
		mary[j].to_int32 = (john[j].to_int32 &  optr) | (mary[j].to_int32 & ~optr);
		john[j].to_int32 = temp;
	}

	if(pnt_b > 0){
		optr=1;
		for(i=0;i<pnt_b-1;i++) optr |= (optr << 1);
		temp = (john[gen_b].to_int32 & ~optr) | (mary[gen_b].to_int32 &  optr);
		mary[gen_b].to_int32 = (john[gen_b].to_int32 &  optr) | (mary[gen_b].to_int32 & ~optr);
		john[gen_b].to_int32 = temp;
	}

	/*
	  printf("john after:\n");
	  print_chrom(john,num_genes,0);
	  printf("mary after:\n");
	  print_chrom(mary,num_genes,0);
	*/

	return;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void mutate(gene *john,int num_genes,double mut_rate){
	/* creates an operator with 1's with rate= mut_rate 
	   uses it to mutate john.
	*/
	int i,j;
	int optr;
	int test;

	for(j=0;j<num_genes;j++){
		optr=0;
		test=1;
		for(i=0;i<32;i++){
			if(RandomDouble() < mut_rate){
				optr |= test;
			}
			test <<= 1;
		}
		john[j].to_int32 ^= optr;
	}
  
	return;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void bin_print(int dec, int len){
	int i,val;
	int test=0;
	int op=1;
	op <<= len-1;
	//printf("op=%u\n",op);
	//printf("dec=%u len=%d\n",dec,len);
	for(i=len-1;i>=0;i--){
		test = (int)pow(2.0f,i);
		//printf("\n[%u]&[%u]=%u test=%u: ",dec,op,dec&op,test);
		val=0;
		if((dec&op) == test) val=1;
		printf("%1d",val);
		op >>= 1;
	}
	return;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void swap_chrom(chromosome *x, chromosome *y){
	chromosome t=*x;*x=*y;*y=t;
}

void quicksort_evalue(chromosome* list,int m,int n){
	double key;
	int i,j;

	//int k;
	if( m < n ) {

		key = list[m].evalue;
		i = m+1;
		j = n;

		while(i <= j)
		{
			while((i <= n) && (list[i].evalue > key))
				i++;
			while((j > m) && (list[j].evalue <= key))
				j--;
			if( i < j)
				swap_chrom(&list[i],&list[j]);
		}

		// swap two elements
		swap_chrom(&list[m],&list[j]);
		// recursively sort the lesser list
		quicksort_evalue(list,m,j-1);
		quicksort_evalue(list,j+1,n);
	}
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
//decreasing order
void quicksort_fitnes(chromosome* list,int m,int n){
	double key;
	int i,j;
  
	//int k;
	if( m < n ) {
		
		key = list[m].fitnes;
		i = m+1;
		j = n;

		while(i <= j)
		{
			while((i <= n) && (list[i].fitnes > key))
				i++;
			while((j > m) && (list[j].fitnes <= key))
				j--;
			if( i < j)
				swap_chrom(&list[i],&list[j]);
		}

		// swap two elements
		swap_chrom(&list[m],&list[j]);
		// recursively sort the lesser list
		quicksort_fitnes(list,m,j-1);
		quicksort_fitnes(list,j+1,n);
	}
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
//decreasing order
int remove_dups(chromosome* list, int num_chrom, int num_genes){
	
	int at=1;
	
	for(int i=1; i<num_chrom; i++){
		//int cmp_chrom2pop_int(const chromosome* chrom,const gene* genes, int num_genes,int start, int last){
		if(cmp_chrom2pop_int(list,list[i].genes,num_genes,0,at)==0){
			swap_chrom(&list[at++],&list[i]);
		}
	}

	return at;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void print_par(const chromosome* chrom,const genlim* gene_lim,int num_chrom,int num_genes){

	for(int i=0;i<num_chrom;i++){
		printf("%2d (",i);
		for(int j=0;j<num_genes;j++) printf("%10.2f ", chrom[i].genes[j].to_ic);
		printf(") ");
		printf(" value=%9.3f fitnes=%9.3f\n",chrom[i].evalue,chrom[i].fitnes);
    
		fflush(stdout);
	}
	return;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void write_par(const chromosome* chrom,const genlim* gene_lim,int ger,char* outfile,int num_chrom,int num_genes){
	int i,j;
	FILE *outfile_ptr;
  
	outfile_ptr=NULL;
	if(!OpenFile_B(outfile,"w",&outfile_ptr)){
		Terminate(6);
	}else{
		fprintf(outfile_ptr,"Generation: %3d\n",ger);

		for(i=0;i<num_chrom;i++){      
			fprintf(outfile_ptr,"%2d (",i);
			for(j=0;j<num_genes;j++) fprintf(outfile_ptr,"%10.2f ",chrom[i].genes[j].to_ic);
			fprintf(outfile_ptr,") ");
			fprintf(outfile_ptr," value=%9.3f fitnes=%9.3f\n",chrom[i].evalue,chrom[i].fitnes);
		}
	}

	CloseFile_B(&outfile_ptr,"w");

	return;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void print_pop(const chromosome* chrom,const genlim* gene_lim,int numc, int numg){
	int i,j;

	for(i=0;i<numc;i++){
		printf("%2d (",i);
		for(j=0;j<numg;j++){printf(" %10d",chrom[i].genes[j].to_int32);}
		printf(") ");
		for(j=0;j<numg;j++){printf(" "),bin_print(chrom[i].genes[j].to_int32,(MAX_GEN_LENGTH));}
		//printf(" value=%10.3f fitnes=%5.3f\n",chrom[i].evalue,chrom[i].fitnes);
		printf("\n");
	}
	return;
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void print_chrom(const chromosome* chrom, int num_genes, int real_flag){
	int j;
	//int i;

	printf("(");
	for(j=0;j<num_genes;j++){
		if(real_flag){
			printf(" %10.5f",chrom->genes[j].to_ic);
		}else{
			printf(" %10d",chrom->genes[j].to_int32);
		}
	}
	printf(") ");
	printf("\n");

	return;
}

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void print_chrom(const gene* genes, int num_genes, int real_flag){
	int j;
	//int i;

	printf("(");
	for(j=0;j<num_genes;j++){
		if(real_flag){
			printf(" %10.5f",genes[j].to_ic);
		}else{
			printf(" %10d",genes[j].to_int32);
		}
	}
	printf(") ");
	printf("\n");

	return;
}

/********************************************************************************
 * This function calculates the RSMD between atomic coordinates of the atoms in *
 * the register ori_ligatm and those for the atoms of the ligand in             *
 * residue[opt_res[0]] after reconstructing the coordinates using opt_par       *
 ********************************************************************************/

double calc_rmsp(int npar, const gene* g1, const gene* g2){
	double rmsp=0.0;
	int i;

	for(i=0;i<npar;i++){
		rmsp += (g1[i].to_ic-g2[i].to_ic)*(g1[i].to_ic-g2[i].to_ic);
	}

	rmsp /= (double)npar;

	return sqrt(rmsp);
}

double genetoic(const genlim* gene_lim, boost::int32_t gene){

	int i=0;
	double tot=gene_lim->bin;
	
	while(tot < RandomDouble(gene)){
		tot += gene_lim->bin;
		i++;
	}

	double ic = gene_lim->min + gene_lim->del * (double)i;
	
	/*        
	printf("ic=%.1f gene=%d randdouble=%.8f min=%.3f del=%.3f bin=%.8f\n", 
	       ic, gene, RandomDouble(gene), 
	       gene_lim->min, gene_lim->del, gene_lim->bin);
	*/

	return(ic);
}

int ictogene(const genlim* gene_lim, double ic){

	int i = (int)((ic - gene_lim->min) / gene_lim->del);

	double tot = 1.0;
	
	while(i > 0){
		tot -= gene_lim->bin;
		i--;
	}

	int gene = RandomInt(tot);
	
        /*
	printf("ic=%.3f gene=%d randdouble=%.5f min=%.3f del=%.3f bin=%.3f\n", 
	       ic, gene, RandomDouble(gene), 
	       gene_lim->min, gene_lim->del, gene_lim->bin);
	*/
	
	return(gene);
}


int RandomInt(double frac){
	return (int)(frac*((double)RAND_MAX+1.0));
}

double RandomDouble(boost::int32_t gene){
	return gene/((double)MAX_RANDOM_VALUE+1.0);
}

double RandomDouble(){
	return rand()/((double)RAND_MAX+1.0);
}
