#include "gaboom.h"
#include "boinc.h"

void cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
	bool Hungarian = false;

	int i,j;
	cfstr cf;                                /* complementarity function value */
	resid *res_ptr = NULL;
	cfstr* cf_ptr = NULL;

	FILE* outfile_ptr = NULL;

	// will 
	double partition_function = 0.0;

	float rmsd = 0.0f;
	int num_of_results = FA->max_results;
	int num_of_clusters = 0;
	int n_unclus = 0;

	char sufix[10];
	char remark[MAX_REMARK];
	char tmpremark[MAX_REMARK];

	// Clustering Variable Definitions
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
	Clus_GAPOP = (int*)malloc(num_chrom*sizeof(int));		// Population (Cluster)	
	Clus_RMSDT = (float*)malloc(num_chrom*sizeof(float));	// RMSD
	Clus_ACF = (double*)malloc(num_chrom*sizeof(double));	// Apparent CF
	Clus_TCF = (double*)malloc(num_chrom*sizeof(double));	// Total CF
	Clus_TOP = (int*)malloc(num_chrom*sizeof(int));			// Best Chromosome in Cluster
	Clus_FRE = (int*)malloc(num_chrom*sizeof(int));			// Frequency
      
	if(!Clus_GAPOP || !Clus_RMSDT || !Clus_ACF ||
	   !Clus_TCF   || !Clus_TOP   || !Clus_FRE)   
	{
		fprintf(stderr,"ERROR: memory allocation error for clusters\n");
		Terminate(2);
	}
      
    ////////////////////////////////
    //////       END         ///////
    ////////////////////////////////
  
    /******************************************************************/
  
    //-------------------------------------------------------
    // fixed center clustering of chrmosome population around highest ranking
    // solutions with rmsd_threshold angstrons threshold
    // Clus_GAPOP[i]=j assigns for each chromosome i to which cluster it belongs
    // as described by j, the chromosome index of the "cluster head", i.e.,that with the
    // higest CF value.
	n_unclus=num_chrom;
	num_of_clusters=0;
	
	// CLustering Variable Initialization and partition_function calculation
	for(j=0;j<num_chrom;++j)
	{
		Clus_GAPOP[j]=-1;
		Clus_ACF[j]=0.0;
		Clus_TCF[j]=0.0;
		Clus_TOP[j]=0;
		Clus_FRE[j]=0;
		if(FA->temperature){
			partition_function += pow( E, ((-1.0) * FA->beta * chrom[j].app_evalue) );
		}
	}
    //printf("n_unclus=%d\n",n_unclus);
    //PAUSE;
	
	// Verify that partition_function != NULL
	if(FA->temperature && partition_function == 0.0) 
	{
		fprintf(stderr,"ERROR: The Partition Function is NULL in the clustering step.\n");
		Terminate(2);
	}
	
	// Clustering part
	while(n_unclus > 0)
	{
		for(j=0;j<num_chrom;++j){if(Clus_GAPOP[j]==-1){break;}}
		//printf("at chromosome j=%d with app_evalue=%.3f\n", j, chrom[j].app_evalue);
		Clus_GAPOP[j]=j;
		Clus_RMSDT[j]=0.0;
		n_unclus--;
		// Clus_TCF[num_of_clusters]=chrom[j].app_evalue;
		// Clus_TCF[num_of_clusters]=chrom[j].app_evalue;
		if(FA->temperature){
			double Pj = pow( E, ((-1.0) * FA->beta * chrom[j].app_evalue) ) / partition_function;
			Clus_ACF[num_of_clusters] = (double)( ( Pj * chrom[j].app_evalue) - (FA->temperature * Pj * log(Pj)) );
			Clus_TCF[num_of_clusters] = (double)( ( Pj * chrom[j].app_evalue) - (FA->temperature * Pj * log(Pj)) );
		}else{
			Clus_TCF[num_of_clusters] = chrom[j].app_evalue;
			Clus_ACF[num_of_clusters] = chrom[j].app_evalue;
		}
		Clus_TOP[num_of_clusters]=j;
		Clus_FRE[num_of_clusters++];

		//printf("n_unclus=%d j=%d\n",n_unclus,j);
		//PAUSE;
		for(i=j+1;i<num_chrom;++i)
		{
			if(Clus_GAPOP[i]==-1)
			{
				rmsd = calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,GB->num_genes,i,j, NULL, NULL, true);
				//printf("rmsd(%d,%d)=%f\n",i,j,rmsd);
				//PAUSE;
				if(rmsd <= FA->cluster_rmsd)
				{
					Clus_GAPOP[i]=j;
					Clus_RMSDT[i]=rmsd;
					// printf("i=%d belongs to cluster of j=%d because rmsd=%.3f\n", i, j, rmsd);
					n_unclus--;
					if(FA->temperature){
						double Pi = pow( E, ((-1.0) * FA->beta * chrom[i].app_evalue) ) / partition_function;
						Clus_ACF[num_of_clusters] += (double)( (Pi * chrom[i].app_evalue) - (FA->temperature * Pi * log(Pi)) );
					}else{
						Clus_ACF[num_of_clusters] += chrom[i].app_evalue;
					}
					Clus_FRE[num_of_clusters++];
				}
			}
		}
		num_of_clusters++;
        
		// quit storing clusters up to N max results
		if(num_of_clusters == num_of_results){break;}
	}

	if(FA->temperature){
		// Reordering the clusters properly by lowest ACF values first (after considering cluster's entropy !)
		QuickSort_Clusters(Clus_TOP, Clus_FRE, Clus_TCF, Clus_ACF, Clus_GAPOP, 0, num_of_results-1);
	}
      
	// print cluster information
	sprintf(sufix,".cad");
	strcpy(tmp_end_strfile,end_strfile);
	strcat(tmp_end_strfile,sufix);
      
	if(!OpenFile_B(tmp_end_strfile,"w",&outfile_ptr))
	{
		Terminate(6);
	}
	else
	{
		for(i=0;i<num_of_clusters;++i)
		{
			fprintf(outfile_ptr,"Cluster %d: TOP=%d TCF=%f ACF=%f freq=%d\n",i,
				Clus_TOP[i],Clus_TCF[i],
				Clus_ACF[i], Clus_FRE[i]);
		}
		if(num_of_clusters > 1)
		{
			fprintf(outfile_ptr,"RMSD between clusters\n");
			for(i=0;i<num_of_clusters;++i)
			{
				for(j=i+1;j<num_of_clusters;++j)
				{
					rmsd=calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,GB->num_genes,Clus_TOP[i],Clus_TOP[j], NULL, NULL, true);
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
      
	for(j=0;j<num_of_results;++j)
	{
		// get parameters of fittest individual in population
		// after optimization -> best docking candidate
    
		// cf=chrom[Clus_TOP[j]].app_evalue;

		for(int k=0; k<GB->num_genes; ++k)
		{
			FA->opt_par[k] = chrom[Clus_TOP[j]].genes[k].to_ic;
		}

		cf=ic2cf(FA,VC,atoms,residue,cleftgrid,GB->num_genes,FA->opt_par);

		strcpy(remark,"REMARK optimized structure\n");

		sprintf(tmpremark,"REMARK CF=%8.5f\n",get_cf_evalue(&cf));
		strcat(remark,tmpremark);
		sprintf(tmpremark,"REMARK CF.app=%8.5f\n",get_apparent_cf_evalue(&cf));
		strcat(remark,tmpremark);

		for(i=0;i<FA->num_optres;++i)
		{
	  
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
			sprintf(tmpremark,"REMARK Residue has an overall SAS of %.3f\n",cf_ptr->totsas);
			strcat(remark,tmpremark);
		}

		sprintf(tmpremark,"REMARK Cluster %d: Rank (top):%d Average CF:%8.5f Frequency:%d\n",
			j,Clus_TOP[j],Clus_ACF[j],Clus_FRE[j]);
		strcat(remark,tmpremark);
		for(i=0;i<FA->npar;++i)
		{
			sprintf(tmpremark,"REMARK [%8.3f]\n",FA->opt_par[i]);
			strcat(remark,tmpremark);
		}
		//sprintf(tmpremark,"REMARK seed=%ld\n",FA->seed_ini);
		strcat(remark,tmpremark);
		if(FA->refstructure == 1){
			sprintf(tmpremark,"REMARK %8.5f RMSD to ref. structure (no symmetry correction)\n",
				calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian));
			strcat(remark,tmpremark);
			Hungarian = true;
			sprintf(tmpremark,"REMARK %8.5f RMSD to ref. structure     (symmetry corrected)\n",
				calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian));
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
	if(FA->refstructure == 1){write_rrd(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,Clus_GAPOP,Clus_RMSDT,end_strfile); }


	// Clusters memory de-allocation
	if(Clus_GAPOP != NULL) free(Clus_GAPOP);
	if(Clus_RMSDT != NULL) free(Clus_RMSDT);
	if(Clus_ACF   != NULL) free(Clus_ACF);
	if(Clus_TCF   != NULL) free(Clus_TCF);
	if(Clus_TOP   != NULL) free(Clus_TOP);
	if(Clus_FRE   != NULL) free(Clus_FRE);
	
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*                  QuickSort functions for Clusters                   */
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void QuickSort_Clusters(int* TOP, int* FRE, double* TCF, double* ACF, int* GAPOP, int beg, int end)
{
	int l,r,p;
	double piv;

	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		piv = ACF[p];
		
		while(1)
		{
			while( (l<=r) && QS_ASC(ACF[l],piv) <= 0 ) ++l;
			while( (l<=r) && QS_ASC(ACF[r],piv)  > 0 ) --r;
			
			if (l > r) break;
			
			swap_clusters(&TOP[l], &FRE[l], &TCF[l], &ACF[l], &GAPOP[l],&TOP[r], &FRE[r], &TCF[r], &ACF[r], &GAPOP[r]);
			
			if (p == r) p=l;
			++l;--r;
		}
		swap_clusters(&TOP[p], &FRE[p], &TCF[p], &ACF[p], &GAPOP[p],&TOP[r], &FRE[r], &TCF[r], &ACF[r], &GAPOP[r]);
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_Clusters(TOP, FRE, TCF, ACF, GAPOP, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_Clusters(TOP, FRE, TCF, ACF, GAPOP, l, end);
			end = r;
		}
	}
}
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*                   Swap Function for Clusters                        */
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
void swap_clusters(int* TOPx, int* FREx, double* TCFx, double* ACFx, int* GAPOPx, int* TOPy, int* FREy, double* TCFy, double* ACFy, int* GAPOPy)
{
	int TOPt, FREt, GAPOPt;
	double TCFt, ACFt;
	TOPt = *TOPx; *TOPx = *TOPy; *TOPy = TOPt;
	FREt = *FREx; *FREx = *FREy; *FREy = FREt;
	TCFt = *TCFx; *TCFx = *TCFy; *TCFy = TCFt;
	ACFt = *ACFx; *ACFx = *ACFy; *ACFy = ACFt;
	GAPOPt = *GAPOPx; *GAPOPx = *GAPOPy; *GAPOPy = GAPOPt;
}
