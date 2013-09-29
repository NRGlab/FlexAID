#include "gaboom.h"
#include "boinc.h"

void cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, 
	     gridpoint* cleftgrid, int num_chrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp){

	int i,j;
	cfstr cf;                                /* complementarity function value */
	resid *res_ptr = NULL;
	cfstr* cf_ptr = NULL;

	FILE* outfile_ptr = NULL;

	float rmsd = 0.0f;
	int num_of_results = FA->max_results;
	int num_of_clusters = 0;
	int n_unclus = 0;

	char sufix[10];
	char remark[MAX_REMARK];
	char tmpremark[MAX_REMARK];

	int* Clus_GAPOP = NULL;
	float* Clus_RMSDT = NULL;
	double* Clus_ACF = NULL;
	double* Clus_TCF = NULL;
	int* Clus_TOP = NULL;
	int* Clus_FRE = NULL;

        ////////////////////////////////
        //// memory allocation for /////
        //// for clustering chrom  /////
        ////////////////////////////////
	Clus_GAPOP = (int*)malloc(num_chrom*sizeof(int));
	Clus_RMSDT = (float*)malloc(num_chrom*sizeof(float));
	Clus_ACF = (double*)malloc(num_chrom*sizeof(double));
	Clus_TCF = (double*)malloc(num_chrom*sizeof(double));
	Clus_TOP = (int*)malloc(num_chrom*sizeof(int));
	Clus_FRE = (int*)malloc(num_chrom*sizeof(int));
      
	if(!Clus_GAPOP || !Clus_RMSDT || !Clus_ACF ||
	   !Clus_TCF   || !Clus_TOP   || !Clus_FRE)   {
		fprintf(stderr,"ERROR: memory allocation error for clusters\n");
		Terminate(2);
	}
      
        ////////////////////////////////
        //////       END         ///////
        ////////////////////////////////
      
        /******************************************************************/
      
        //---------------------------------------------------------
        // fixed center clustering of chrmosome population around highest ranking
        // solutions with rmsd_threshold angstrons threshold
        // Clus_GAPOP[i]=j assigns for each chromosome i to which cluster it belongs
        // as described by j, the chromosome index of the "cluster head", i.e.,that with the
        // higest CF value.
	n_unclus=num_chrom;
	num_of_clusters=0;
	for(j=0;j<num_chrom;j++){
		Clus_GAPOP[j]=-1;
		Clus_ACF[j]=0.0;
		Clus_TCF[j]=0.0;
		Clus_TOP[j]=0;
		Clus_FRE[j]=0;
	}
        //printf("n_unclus=%d\n",n_unclus);
        //PAUSE;
	
	while(n_unclus > 0){
		for(j=0;j<num_chrom;j++){if(Clus_GAPOP[j]==-1){break;}}
		//printf("at chromosome j=%d with app_evalue=%.3f\n", j, chrom[j].app_evalue);
		Clus_GAPOP[j]=j;
		Clus_RMSDT[j]=0.0;
		n_unclus--;
		Clus_TCF[num_of_clusters]=chrom[j].app_evalue;
		Clus_ACF[num_of_clusters]=chrom[j].app_evalue;
		Clus_TOP[num_of_clusters]=j;
		Clus_FRE[num_of_clusters]++;
		//printf("n_unclus=%d j=%d\n",n_unclus,j);
		//PAUSE;
		for(i=j+1;i<num_chrom;i++){
			if(Clus_GAPOP[i]==-1){
				rmsd=calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,GB->num_genes,i,j);
				//printf("rmsd(%d,%d)=%f\n",i,j,rmsd);
				//PAUSE;
				if(rmsd<=FA->cluster_rmsd){
					Clus_GAPOP[i]=j;
					Clus_RMSDT[i]=rmsd;
					//printf("i=%d belongs to cluster of j=%d because rmsd=%.3f\n", i, j, rmsd);
					n_unclus--;
					Clus_ACF[num_of_clusters] += chrom[i].app_evalue;
					Clus_FRE[num_of_clusters]++;
				}
			}
		}
		Clus_ACF[num_of_clusters] /= (double)Clus_FRE[num_of_clusters];
		num_of_clusters++;
        
		// quit storing clusters up to N max results
		if(num_of_clusters == num_of_results){break;}
	}
      
	// print cluster information
	sprintf(sufix,".cad");
	strcpy(tmp_end_strfile,end_strfile);
	strcat(tmp_end_strfile,sufix);
      
	if(!OpenFile_B(tmp_end_strfile,"w",&outfile_ptr)){
		Terminate(6);
	}else{
		for(i=0;i<num_of_clusters;i++){
			fprintf(outfile_ptr,"Cluster %d: TOP=%d TCF=%f ACF=%f freq=%d\n",i,
				Clus_TOP[i],Clus_TCF[i],
				Clus_ACF[i], Clus_FRE[i]);
		}
		if(num_of_clusters > 1){
			fprintf(outfile_ptr,"RMSD between clusters\n");
			for(i=0;i<num_of_clusters;i++){
				for(j=i+1;j<num_of_clusters;j++){
					rmsd=calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,GB->num_genes,Clus_TOP[i],Clus_TOP[j]);
					fprintf(outfile_ptr,"rmsd(%d,%d)=%f\n",i,j,rmsd);
				}
			} 
		}
	}
	CloseFile_B(&outfile_ptr,"w");
	

	if(num_of_clusters < num_of_results){num_of_results=num_of_clusters;}
        //num_of_results=1;
      
	printf("num_of_clusters=%d num_of_results=%d\n",num_of_clusters,num_of_results);
	
        // output results, 10% of the number of chromosomes or 
        // the number of clusters, the smallest.
      
	for(j=0;j<num_of_results;j++){
		// get parameters of fittest individual in population
		// after optimization -> best docking candidate
    
		// cf=chrom[Clus_TOP[j]].app_evalue;

		for(int k=0; k<GB->num_genes; k++){
			FA->opt_par[k] = chrom[Clus_TOP[j]].genes[k].to_ic;
		}

		cf=ic2cf(FA,VC,atoms,residue,cleftgrid,GB->num_genes,FA->opt_par);

		strcpy(remark,"REMARK optimized structure\n");

		sprintf(tmpremark,"REMARK CF=%8.5f\n",get_cf_evalue(&cf));
		strcat(remark,tmpremark);
		sprintf(tmpremark,"REMARK CF.app=%8.5f\n",get_apparent_cf_evalue(&cf));
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

		sprintf(tmpremark,"REMARK Cluster %d: Rank (top):%d Average CF:%8.5f Frequency:%d\n",
			j,Clus_TOP[j],Clus_ACF[j],Clus_FRE[j]);
		strcat(remark,tmpremark);
		for(i=0;i<FA->npar;i++){
			sprintf(tmpremark,"REMARK [%8.3f]\n",FA->opt_par[i]);
			strcat(remark,tmpremark);
		}
		//sprintf(tmpremark,"REMARK seed=%ld\n",FA->seed_ini);
		strcat(remark,tmpremark);
		if(FA->refstructure == 1){
			sprintf(tmpremark,"REMARK %8.5f RMSD to ref. structure\n",
				calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par));
			strcat(remark,tmpremark);
		}
		sprintf(tmpremark,"REMARK inputs: %s & %s\n",dockinp,gainp);
		strcat(remark,tmpremark);
		sprintf(sufix,"_%d.pdb",j);
		strcpy(tmp_end_strfile,end_strfile);
		strcat(tmp_end_strfile,sufix);
		//printf("filename=<%s>\n",tmp_end_strfile);
		//PAUSE;
		write_pdb(FA,atoms,residue,tmp_end_strfile,remark);
	}
      
        // print the RMSD between each chrom. and the reference structure if there is one.
	if(FA->refstructure == 1){ write_rrd(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,Clus_GAPOP,Clus_RMSDT,end_strfile); }


	// Clusters memory de-allocation
	if(Clus_GAPOP != NULL) free(Clus_GAPOP);
	if(Clus_RMSDT != NULL) free(Clus_RMSDT);
	if(Clus_ACF != NULL) free(Clus_ACF);
	if(Clus_TCF != NULL) free(Clus_TCF);
	if(Clus_TOP != NULL) free(Clus_TOP);
	if(Clus_FRE != NULL) free(Clus_FRE);
	
}
