#include "gaboom.h"
#include "boinc.h"

// Chi(x) function as defined in :
//   Science 344(6191):1492-1496, 2014, Eq. (1)
#define Chi(a,d) 	( ((a-d) < 0.0) ? 1 : 0 )

// #DEFINE neighborRateLow and neighborRateHigh 
#define NEIGHBORRATELOW 0.01
#define NEIGHBORRATEHIGH 0.02
#define EXCLUDE_HALO true
#define OUTPUT_CLUSTER_CENTER true

void DensityPeak_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
	// Density Peak Clustering variables declaration
	int i,j,k;
	bool Hungarian = false;
	bool Entropic = ( FA->temperature > 0 ? true : false );
	float dc = 0.0f;
	const int nAtoms = residue[atoms[FA->map_par[0].atm].ofres].latm[0] - residue[atoms[FA->map_par[0].atm].ofres].fatm[0] + 1;
	uint maxDensity;
	int mean, stddev;
	int nResults;
	int nClusters, nUnclustered, nOutliers;
	float maxDist, minDist;
	float* RMSD;
	double Pi;
	double partition_function;
	ClusterChrom* Chrom;
	ClusterChrom *pChrom, *iChrom, *jChrom;
	Cluster* Clust;
	Cluster *pCluster, *iClust, *jClust;

	// File and Output variables declarations
	cfstr CF;                                /* complementarity function value */
	resid *pRes = NULL;
	cfstr* pCF = NULL;

	FILE* outfile_ptr = NULL;
	char sufix[10];
	char remark[MAX_REMARK];
	char tmpremark[MAX_REMARK];

	//  dynamically allocated memory check-up
	Chrom = (ClusterChrom*) malloc(num_chrom*sizeof(ClusterChrom));
	if(Chrom == NULL)
	{
		fprintf(stderr,"ERROR: memory allocation error for ChromClusters data structures.\n");
		Terminate(2);
	}

	RMSD = (float*) malloc(num_chrom*num_chrom*sizeof(float));
	if(RMSD == NULL)
	{
		fprintf(stderr,"ERROR: memory allocation error for RMSD matrix.\n");
		Terminate(2);
	}

	// variables initialization
	memset(RMSD, 0.0, num_chrom * num_chrom);
	for(i = 0, partition_function = 0.0; i < num_chrom; ++i)
	{
		pChrom = &Chrom[i];
		if(pChrom)
		{
			pChrom->index = i;
			pChrom->isHalo = false;
			pChrom->isCenter = false;
			pChrom->isBorder = false;
			pChrom->isClustered = false;
			pChrom->Chromosome = &chrom[i];
			pChrom->Cluster = 0;
			pChrom->Density = 0;
			pChrom->CF = 0.0;
			pChrom->DP = NULL;
			pChrom->DPdist = 0.0;
			memset(pChrom->Coord, 0.0, 3*MAX_ATM_HET);
			if(Entropic) partition_function += pow( E, ((-1.0) * FA->beta * pChrom->Chromosome->app_evalue) );
		}
	}

	// (0) VERIFICATIONS 
	// Verify that partition_function != NULL
	if( Entropic && partition_function == 0.0 ) 
	{
		fprintf(stderr,"ERROR: The Partition Function is NULL during the clustering step.\n");
		Terminate(2);
	}

	// (1) Build Chromosome Cartesian Coordinates
	for(i = 0; i < num_chrom; ++i)
	{
		pChrom = &Chrom[i];
		if(i+1 == num_chrom) calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,GB->num_genes,i,i, pChrom->Coord, NULL, false);
        else calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,GB->num_genes,i,i+1, pChrom->Coord, (pChrom++)->Coord, false);
	}

	// (2) Build RMSD Matrix
	for(i = 0, Pi = 0.0, iChrom=NULL; i < num_chrom; ++i)
	{
		iChrom = &Chrom[i];
		if(Entropic)
		{
			Pi = pow( E, ((-1.0) * FA->beta * iChrom->Chromosome->app_evalue) ) / partition_function;
			iChrom->CF = (double) ( Pi * iChrom->Chromosome->app_evalue) - (FA->temperature * Pi * log(Pi));
		}
		else iChrom->CF = iChrom->Chromosome->app_evalue;
		for(j = i+1; j < num_chrom; ++j)
		{
			jChrom = &Chrom[j];
			for(k = 0, minDist=0.0; k < nAtoms; ++k) minDist += sqrdist(&iChrom->Coord[3*k], &jChrom->Coord[3*k]);
			RMSD[K(i,j,num_chrom)] = sqrtf(minDist/(float)nAtoms);
		}
	}

	// (*) Determine Distance Cut-Off
	dc = getDistanceCutoff(RMSD, NEIGHBORRATELOW, NEIGHBORRATEHIGH, num_chrom);
	// dc = 2.5;
	printf("DC:%g\n",dc);

	// (3) Build Local Density Matrix
	for(i = 0; i < num_chrom; ++i)
	{
		iChrom = &Chrom[i];
		// for(j = i+1; j < num_chrom; ++j)
		for(j = 0; j < num_chrom; ++j)
		{
			jChrom	= &Chrom[j];
			if(jChrom != iChrom) iChrom->Density += Chi(RMSD[K(i,j,num_chrom)], dc);
		}
	}

	// (4) Fill out DP and DPdist in Chrom
	for(i = 0; i < num_chrom; ++i)
	{
		iChrom = &Chrom[i];
		for(j = 0, minDist=FLT_MAX; j < num_chrom; ++j)
		{
			if(j != i)
			{
				jChrom = &Chrom[j];
				if(jChrom->Density > iChrom->Density && RMSD[K(i,j,num_chrom)] < minDist && RMSD[K(i,j,num_chrom)] > 0.0)
				{
					minDist = RMSD[K(i,j,num_chrom)];
					iChrom->DP = jChrom;
					iChrom->DPdist = minDist;
				}
			}
		}
		for(j = num_chrom-1, minDist=FLT_MAX; j >= 0; --j)
		{
			if(j != i)
			{
				jChrom = &Chrom[j];
				if(jChrom->Density > iChrom->Density && RMSD[K(i,j,num_chrom)] < minDist && RMSD[K(i,j,num_chrom)] > 0.0)
				{
					minDist = RMSD[K(i,j,num_chrom)];
					iChrom->DP = jChrom;
					iChrom->DPdist = minDist;
				}
			}
		}
	}

	// (5) Find Maximal DPdist value in data (maxDPdist -> maxDist, highest density chrom -> pChrom)
	for(pChrom=NULL, i=0, maxDist=0.0, maxDensity = 0, k=0; i < num_chrom; ++i)
	{
		iChrom = &Chrom[i];
		if(iChrom->DPdist > maxDist) maxDist = iChrom->DPdist;
		if(iChrom->Density > maxDensity)
		{
			pChrom = iChrom;
			maxDensity = pChrom->Density;
		}
		if(iChrom->DP == NULL) k++;
		iChrom->PiDi = iChrom->Density * iChrom->DPdist;
	}

	// (5.1) pChrom (which points to the highest density chrom) should not have pChrom->DP defined
	if(pChrom != NULL)
	{
		if(pChrom->DP == NULL)
		{
			pChrom->DPdist = maxDist;
			pChrom->PiDi = pChrom->Density * pChrom->DPdist;
			--k;
		}
		for(i=0, iChrom=Chrom; i<num_chrom; ++iChrom, ++i)
		{
			if(iChrom == pChrom) continue;
			if(iChrom->DP == NULL && iChrom->Density == maxDensity) 
			{
				if(RMSD[K(iChrom->index,pChrom->index,num_chrom)] <= dc)
				{
					iChrom->DP = pChrom;
					iChrom->DPdist = RMSD[K(iChrom->index,pChrom->index,num_chrom)];
					iChrom->PiDi = iChrom->Density * iChrom->DPdist;
					--k;
				}
				else
				{
					iChrom->DPdist = maxDist;
					iChrom->PiDi = iChrom->Density * iChrom->DPdist;
					--k;
				}
			}
		}
	}
	// (5.2) Dealing with multiple Chrom->DP == NULL
	while(k)
	{
		i = 0; pChrom = Chrom; // equivalent to pChrom = &Chrom[0]
		while(pChrom->DP != NULL /* && pChrom->DPdist > 0.0*/ && i < num_chrom) { ++pChrom; i++; }
		if(pChrom->Density == maxDensity) pChrom->DPdist = maxDist;
		else 
		{
			minDist = FLT_MAX;
			for(j=0, jChrom=NULL; j < num_chrom; ++j)
			{
				jChrom = &Chrom[j];
				if(jChrom != pChrom && jChrom->Density > pChrom->Density && RMSD[K(i,j,num_chrom)] <= minDist && RMSD[K(i,j,num_chrom)] > 0.0)
				{
					minDist = RMSD[K(i,j,num_chrom)];
					pChrom->DP = jChrom;
					pChrom->DPdist = minDist;
				}
			}
		}
		pChrom->PiDi = pChrom->Density * pChrom->DPdist;
		--k;
	}
	

	// (5.3) Sort Chrom by decreasing PiDi value
	QuickSort_ChromCluster_by_PiDi(Chrom,num_chrom,0,num_chrom-1);

	// (6) Identify Cluster Centers
	for(pChrom=NULL, i=0, nClusters=0, stddev=calculate_mean(Chrom, num_chrom), mean=calculate_mean(Chrom, num_chrom); (&Chrom[nClusters])->PiDi > (mean + 2*stddev) && i < num_chrom; ++i)
	{
		pChrom = &Chrom[nClusters];
		
		if(pChrom != NULL)
		{	
			if(pChrom->DP != NULL && pChrom->DP->Cluster >= 1) pChrom->Cluster = pChrom->DP->Cluster;
			else
			{
				pChrom->Cluster = (++nClusters);		// cluster assigned
				pChrom->isCenter = true; 				// cluster center assigned
			}
		}
	}
	printf("nClusters:%d\n",nClusters);
	// (6.1) QuickSort by decreasing Density value
	QuickSort_ChromCluster_by_Density(Chrom, num_chrom, 0, num_chrom-1);

	// (7) Clustering Step
	for(i=0, pChrom=NULL; i<num_chrom; ++i)
	{
		pChrom = &Chrom[i];
		iChrom = &Chrom[i];
		while(pChrom->Cluster <= 0 && pChrom->DP != NULL) pChrom = pChrom->DP;
		if(pChrom->Cluster > 0) iChrom->Cluster = pChrom->Cluster;
	}

	// Sorting ChromCluster elements by ASCENCING CF values
	QuickSort_ChromCluster_by_CF(Chrom, num_chrom, 0, num_chrom-1);

	// (8) Assignation of chromosome to its cluster Core/Halo 
	// At this point, the cluster core vs cluster halo assignation is unuseful if there is a single cluster (unable to separe halore (noise) from core data points in dataset)
	if(EXCLUDE_HALO && nClusters > 1)
	{
		// (8.1) Find for each Cluster(k): define the border region 
		// (8.2) Find for each Cluster(m): the point of highest density(Pb) within the border region
		for(k = 1, maxDensity=0, pChrom=NULL; k <= nClusters; ++k)
		{
			for(i=0, iChrom=Chrom; i<num_chrom; ++i, ++iChrom) if( k == iChrom->Cluster )
			{
				for(j=0, jChrom=Chrom; j<num_chrom; ++j, ++jChrom) if(jChrom->Cluster > 0 && jChrom->Cluster != k)
				{
					if( RMSD[K(iChrom->index, jChrom->index, num_chrom)] < ( (dc < FA->cluster_rmsd) ? dc*0.5 : FA->cluster_rmsd) /*&& iChrom->Density > maxDensity*/)
					{
						iChrom->isBorder = true;
					}
				}
			}
			for(i=0, iChrom = Chrom; i<num_chrom; ++i, ++iChrom) if( k == iChrom->Cluster)
			{
				if(iChrom->isBorder && iChrom->Density > maxDensity)
				{
					pChrom = iChrom;
					maxDensity = pChrom->Density;
				}
			}
            if(pChrom)
            {
            	for(i=0, iChrom=Chrom; i<num_chrom; ++i, ++iChrom) if( k == iChrom->Cluster) if(iChrom->Density < pChrom->Density) iChrom->isHalo = true;
			}
		}
		for(i=0, iChrom=Chrom,nOutliers=0,j=0;i<num_chrom;++i,iChrom++)
		{
			if( iChrom->isHalo && iChrom->Chromosome->app_evalue > 1000) nOutliers++;
			if( iChrom->isHalo ) j++;
		}
		printf("nOutliers:%d\tnHalo:%d\n",nOutliers, j);
	}

	// (9) Cluster Creation
	if(nClusters < FA->max_results) nResults = nClusters;
	else nResults = FA->max_results;

	Clust = (Cluster*) malloc( nResults * sizeof(Cluster) );
	// //  dynamically allocated memory check-up
	if(Clust == NULL) 
	{
		fprintf(stderr,"ERROR: memory allocation error for Clusters data structures\n");
		Terminate(2);
	}
	for(i=0; i < nResults; ++i)
	{
		pCluster = &Clust[i];
		// initializing the Cluster element
		pCluster->ID = 0;
		pCluster->Frequency = 0;
		pCluster->totCF = 0.0;
		pCluster->lowestCF = 0.0;
		pCluster->BestCF = NULL;
		pCluster->Center = NULL;
	}
	// (10) Building UP the Clusters 
	for(k = 1, pCluster=NULL; k <= nResults; ++k)
	{
		pCluster = &Clust[k-1];
		pCluster->ID = k;
		for(i=0; i<num_chrom; ++i)
		{
            pChrom = &Chrom[i];
            if(pChrom && !pChrom->isHalo && pChrom->Cluster == k && !pChrom->isClustered )
            {
                pCluster->totCF += pChrom->CF;
                pCluster->Frequency += 1;
                if( pChrom->isCenter  ) pCluster->Center = pChrom;
                if( !pCluster->BestCF ) { pCluster->BestCF = pChrom; pCluster->lowestCF = pChrom->CF; }
                else if(pChrom->CF < pCluster->BestCF->CF) { pCluster->BestCF = pChrom; pCluster->lowestCF = pChrom->CF; }
                pChrom->isClustered = true;
            }
		}
	}
	
	// looking for an interesting individual chrom to output 
	for(i=0; i < num_chrom && k <= nResults; ++i)
	{
        pChrom = &Chrom[i];
        if (pChrom && !pChrom->isClustered && pChrom->Cluster == 0)// && pChrom->Chromosome->app_evalue < 1000)
		{
			pCluster = &Clust[k];
			pCluster->ID = k;
			pCluster->totCF = pChrom->CF;
			pCluster->Frequency = 1;
			pCluster->Center = pChrom;
			pCluster->BestCF = pChrom;
			pCluster->lowestCF = pChrom->CF;

			pChrom->Cluster = k;
			pChrom->isCenter = true;
            pChrom->isClustered = true;
			++k;
		}
	}

	// (11) Sort Clusters by ascending CF (lowest first) && Chromosomes by DESCENDING Density
    QuickSort_Cluster_by_CF( Clust, Entropic, 0, nResults-1 );

	// (12) Output Cluster information
	sprintf(sufix,".cad");
	strcpy(tmp_end_strfile, end_strfile);
	strcat(tmp_end_strfile, sufix);

	if(!OpenFile_B(tmp_end_strfile,"w",&outfile_ptr))
	{
		Terminate(6);
	}
	else
	{
		for(i = 0; i < nResults; ++i)
		{
            pCluster = &Clust[i];
			if(pCluster && pCluster->BestCF && pCluster->Center) fprintf(outfile_ptr, "Cluster %d: Center:%d (CF:%g) Best:%d (CF:%g) TCF=%f Frequency=%d\n",
				pCluster->ID,
				pCluster->Center->index, 
				pCluster->Center->CF,
				pCluster->BestCF->index,
				pCluster->BestCF->CF,
				pCluster->totCF,
				pCluster->Frequency);
		}

		if(nClusters > 1)
		{
			fprintf(outfile_ptr,"RMSD between clusters\n");
			for(i=0, iClust = NULL; i < nResults; ++i)
			{
				iClust = &Clust[i];
				for(j=i+1; j < nResults; ++j)
				{
					jClust = &Clust[j];
                    if(!Clust[j].BestCF || (OUTPUT_CLUSTER_CENTER && !Clust[j].Center) ) continue;
					if(OUTPUT_CLUSTER_CENTER==true) fprintf(outfile_ptr,"rmsd(%d,%d)=%f\n",i+1,j+1,RMSD[K(iClust->Center->index, jClust->Center->index, num_chrom)]);
					else fprintf(outfile_ptr,"rmsd(%d,%d)=%f\n",i+1,j+1,RMSD[K(iClust->BestCF->index, jClust->BestCF->index, num_chrom)]);
				}
			} 
		}
	}
	CloseFile_B(&outfile_ptr,"w");

	// (13) Output Clusters
	// output clusters informations (looping through each cluster)
	for(i=0, pCluster=NULL, pChrom=NULL; i < nResults; i++)
	{
		pCluster = &Clust[i];
		// setting pChrom to either BestCF or ClusterCenter
		if(OUTPUT_CLUSTER_CENTER) 	pChrom = pCluster->Center;
		else 						pChrom = pCluster->BestCF;

		if(!pChrom)
		{
			fprintf(stderr,"ERROR: pChrom is undefined in results output point (13)\n");
			Terminate(2);
		}
		// outputting cluster center
		for(k=0; k<GB->num_genes ; ++k) FA->opt_par[k] = pChrom->Chromosome->genes[k].to_ic;
		
		CF = ic2cf(FA, VC, atoms, residue, cleftgrid, GB->num_genes, FA->opt_par);
		
		strcpy(remark,"REMARK optimized structure\n");
		sprintf(tmpremark,"REMARK Density Peak clustering algorithm used to output %s as cluster representatives\n", (OUTPUT_CLUSTER_CENTER == true ? "the lowest CF" : "the center of highest density"));
		strcat(remark,tmpremark);
		
		sprintf(tmpremark,"REMARK CF=%8.5f\n",get_cf_evalue(&CF));
		strcat(remark,tmpremark);
		sprintf(tmpremark,"REMARK CF.app=%8.5f\n",get_apparent_cf_evalue(&CF));
		strcat(remark,tmpremark);
		for(j = 0; j < FA->num_optres; ++j)
		{
			pRes = &residue[FA->optres[j].rnum];
			pCF  = &FA->optres[j].cf;

			sprintf(tmpremark,"REMARK optimizable residue %s %c %d\n", pRes->name, pRes->chn, pRes->number);
			strcat(remark,tmpremark);
	  
			sprintf(tmpremark ,"REMARK CF.com=%8.5f\n", pCF->com);
			strcat(remark, tmpremark);
			sprintf(tmpremark ,"REMARK CF.sas=%8.5f\n", pCF->sas);
			strcat(remark, tmpremark);
			sprintf(tmpremark ,"REMARK CF.wal=%8.5f\n", pCF->wal);
			strcat(remark, tmpremark);
			sprintf(tmpremark ,"REMARK CF.con=%8.5f\n", pCF->con);
			strcat(remark, tmpremark);
			sprintf(tmpremark, "REMARK Residue has an overall SAS of %.3f\n", pCF->totsas);
			strcat(remark, tmpremark);
		}
		
		sprintf(tmpremark,"REMARK Cluster:%d Best CF in Cluster:%8.5f Cluster Center (CF):%8.5f Cluster Total CF:%8.5f Cluster Frequency:%d\n", 
				pCluster->ID, pCluster->BestCF->CF, pCluster->Center->CF, pCluster->totCF, pCluster->Frequency);
		strcat(remark,tmpremark);
		
		for(j=0; j < FA->npar; ++j)
		{
			sprintf(tmpremark, "REMARK [%8.3f]\n",FA->opt_par[j]);
			strcat(remark,tmpremark);
		}

		// Calculate RMSD value to REFERENCE pose IF REFERENCE is defined in input 
		if(FA->refstructure == 1)
		{
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
		sprintf(sufix,"_%d.pdb",i);
		strcpy(tmp_end_strfile,end_strfile);
		strcat(tmp_end_strfile,sufix);

		// (*) write pdb file
		write_pdb(FA,atoms,residue,tmp_end_strfile,remark);
	}

	// Need to modify write_rrd.c OR  
	if(FA->refstructure == 1) { write_DensityPeak_rrd(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,Chrom,Clust,RMSD,end_strfile); }

	// (*) Printing ChromCluster informations
    for(i=0,j=0;i<num_chrom;++i)
    {
        printf("i:%d\tAdd:%p\tDensity:%d\tDistance:%g\tCluster:%d\tDP:%p\tPiDi:%g\tCF:%g\tisCore:%s\tisCenter:%s\n",(&Chrom[i])->index, &Chrom[i], (&Chrom[i])->Density, (&Chrom[i])->DPdist, (&Chrom[i])->Cluster, (&Chrom[i])->DP, (&Chrom[i])->PiDi, (&Chrom[i])->CF, ((&Chrom[i])->isHalo ? "false" : "true"), ((&Chrom[i])->isCenter ? "true" : "false"));
    }
	
	// (*) Memory deallocation
	if(Chrom != NULL) { free(Chrom); Chrom=NULL; } 
	if(RMSD  != NULL) { free(RMSD);  RMSD=NULL;  }
    if(Clust != NULL) { free(Clust); Clust=NULL; }
}
float getDistanceCutoff(float* RMSD, float neighborRateLow, float neighborRateHigh, int num_chrom)
{
	int i,j;
	float dc = 0.0f;
	int neighbors = 0;
	int nLow = neighborRateLow * num_chrom * num_chrom;
	int nHigh = neighborRateHigh * num_chrom * num_chrom;
	while(neighbors < nLow || neighbors > nHigh)
	{
		neighbors = 0;
		for(i=0; i<num_chrom-1; ++i)
		{
			for(j=0; j<num_chrom; ++j)
			{
				if(i==j) continue;
				if(RMSD[K(i,j,num_chrom)] <= dc) ++neighbors;
				// if(neighbors > nHigh) dc += 0.25;
				// if(neighbors > nHigh) goto DCPLUS;
				if(neighbors > nHigh) goto DCPLUS;
			}
		}
		DCPLUS: dc += 0.01;
	}
	return dc;
}
void QuickSort_Cluster_by_CF(Cluster* Clust, bool Entropic, int beg, int end)
{
	int l, r, p;
	double pivot;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		if(Entropic) pivot = (&Clust[p])->totCF;
        else 		 pivot = (&Clust[p])->lowestCF;

		while(1)
		{
			if(Entropic)
			{
				while( (l<=r) && QS_ASC( (&Clust[l])->totCF, pivot ) <= 0.0 ) ++l;
				while( (l<=r) && QS_ASC( (&Clust[r])->totCF, pivot )  > 0.0 ) --r;
			}
			else
			{
				while( (l<=r) && QS_ASC( (&Clust[l])->lowestCF, pivot ) <= 0.0 ) ++l;
				while( (l<=r) && QS_ASC( (&Clust[r])->lowestCF, pivot )  > 0.0 ) --r;
			}
			if (l > r) break;

			swap_clusters(&Clust[l],&Clust[r]);
			if (p == r) p=l;
			++l;--r;
		}
		swap_clusters(&Clust[p], &Clust[r]);
		
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_Cluster_by_CF(Clust, Entropic, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_Cluster_by_CF(Clust, Entropic, l, end);
			end = r;
		}
	}
}
void swap_clusters(Cluster* xCluster, Cluster* yCluster)
{
	Cluster tCluster = *xCluster; *xCluster = *yCluster; *yCluster = tCluster;
}
void QuickSort_ChromCluster_by_CF(ClusterChrom* Chrom, int num_chrom, int beg, int end)
{
	int l, r, p;
	double pivot;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		pivot = (&Chrom[p])->CF;

		while(1)
		{
			while( (l<=r) && QS_ASC((&Chrom[l])->CF,pivot) <= 0.0 ) ++l;
			while( (l<=r) && QS_ASC((&Chrom[r])->CF,pivot)  > 0.0 ) --r;
			
			if (l > r) break;

			swap_elements(Chrom, &Chrom[l],&Chrom[r], num_chrom);
			if (p == r) p=l;
			++l;--r;
		}
		swap_elements(Chrom, &Chrom[p], &Chrom[r], num_chrom);
		
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_ChromCluster_by_CF(Chrom, num_chrom, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_ChromCluster_by_CF(Chrom, num_chrom, l, end);
			end = r;
		}
	}
}
void QuickSort_ChromCluster_by_Density(ClusterChrom* Chrom, int num_chrom, int beg, int end)
{
	int l, r, p;
	int pivot;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		pivot = (&Chrom[p])->Density;

		while(1)
		{
			while( (l<=r) && QS_DSC((&Chrom[l])->Density,pivot) <= 0 ) ++l;
			while( (l<=r) && QS_DSC((&Chrom[r])->Density,pivot)  > 0 ) --r;
			
			if (l > r) break;

			swap_elements(Chrom, &Chrom[l],&Chrom[r], num_chrom);
			if (p == r) p=l;
			++l;--r;
		}
		swap_elements(Chrom, &Chrom[p], &Chrom[r], num_chrom);
		
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_ChromCluster_by_Density(Chrom, num_chrom, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_ChromCluster_by_Density(Chrom, num_chrom, l, end);
			end = r;
		}
	}
}
void QuickSort_ChromCluster_by_PiDi(ClusterChrom* Chrom, int num_chrom, int beg, int end)
{
	int l, r, p;
	float pivot;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		pivot = (&Chrom[p])->PiDi;

		while(1)
		{
			while( (l<=r) && QS_DSC((&Chrom[l])->PiDi,pivot) <= 0 ) ++l;
			while( (l<=r) && QS_DSC((&Chrom[r])->PiDi,pivot)  > 0 ) --r;
			
			if (l > r) break;

			swap_elements(Chrom, &Chrom[l],&Chrom[r], num_chrom);
			if (p == r) p=l;
			++l;--r;
		}
		swap_elements(Chrom, &Chrom[p], &Chrom[r], num_chrom);
		
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_ChromCluster_by_PiDi(Chrom, num_chrom, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_ChromCluster_by_PiDi(Chrom, num_chrom, l, end);
			end = r;
		}
	}
}
void swap_elements(ClusterChrom* Chrom, ClusterChrom* ChromX, ClusterChrom* ChromY, int num_chrom)
{
	int i;
	ClusterChrom* pChrom;
	ClusterChrom ChromT = *ChromX;*ChromX = *ChromY; *ChromY = ChromT;
	
	for(i=0, pChrom=NULL;i<num_chrom;++i)
	{
		pChrom = &Chrom[i];
		if(pChrom->DP == ChromX) pChrom->DP = ChromY;
		else if(pChrom->DP == ChromY) pChrom->DP = ChromX;
	}
		
}
float calculate_stddev(ClusterChrom* Chrom, int num_chrom)
{
	int i;
	float mean = calculate_mean(Chrom, num_chrom);
	float sqrtot = 0.0f;
    for(i = 0; i < num_chrom; ++i) sqrtot += pow( (mean - (&Chrom[i])->PiDi),2.0 );
	return sqrtf(sqrtot/(float)num_chrom);
}
float calculate_mean(ClusterChrom* Chrom, int num_chrom)
{
	int i;
	float mean = 0.0f;
    for(i = 0; i < num_chrom; ++i) mean += (&Chrom[i])->PiDi;
	return mean /= (float)num_chrom;
}