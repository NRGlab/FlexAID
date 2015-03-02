#include "gaboom.h"
#include "boinc.h"

// Chi(x) function as defined in :
//   Science 344(6191):1492-1496, 2014, Eq. (1)
#define Chi(a,d) 	( ((a-d) < 0.0) ? 1 : 0 )
#define K(i,j,n) ( (i < j) ? (i*n+j) : (j*n+i) )

struct ClusterChrom
{
	uint index;						// original index in chrommose* chrom (at the input of density_cluster function)
	bool isCore;					// is part of cluster core?
	bool isCenter;					// is cluster center?
	chromosome* Chromosome;			// Chromosomes list
	int Cluster;					// Assigned Cluster
	uint Density;					// Density of points in distance cut-off
	double CF;						// Complementarity Function value
	float PiDi;						// Density x DPdist
	float DPdist;					// Nearest highest density peak distance
	float Coord[3*MAX_ATM_HET];		// Cartesian Coordinates
	struct ClusterChrom* DP;		// Nearest Density Peak (point of higher density)
}; typedef struct ClusterChrom ClusterChrom;

// struct ClusterChromQ				// serves as a ClusterChrom Queue data structure
// {
// 	ClusterChrom* cChrom;			// pointer to ClusterChrom data structure
// 	ClusterChrom* Next;				// pointer to next ClusterChrom in Cluster
// }; typedef struct ClusterChromQ;

struct Cluster
{
       int ID;									// assigned cluster number (ID)
       int Frequency;							// observation frequency of this cluster (number of representatives in cluster)
       float refRMSD;							// Cluster RMSD to reference value (if a reference in defined)
       double totCF;				
       // Pointer to best CF value in cluster
       ClusterChrom* BestCF;					// Pointer to the ClusterChrom individual with the lowest CF in cluster
       // Pointer to Cluster Center
       ClusterChrom* Center;					// Queue of ClusterChrom (first element is the cluster center)
};
typedef struct Cluster Cluster;

// struct ClusteredChromosomes
// {
// 	int nClusters;
// 	int nOutliers;
// 	Cluster* first;
// 	Cluster* next;
// 	ClusterChrom* outliers;
// }; typedef struct ClusteredChromosomes ClusteredChrom;

void QuickSort_ChromCluster_by_Density(ClusterChrom* Chrom, int num_chrom, int beg, int end);
void QuickSort_ChromCluster_by_PiDi(ClusterChrom* Chrom, int num_chrom, int beg, int end);
void swap_elements(ClusterChrom* Chrom, ClusterChrom* ChromX, ClusterChrom* ChromY, int num_chrom);
float calculate_stddev(ClusterChrom* Chrom, int num_chrom);
float calculate_mean(ClusterChrom* Chrom, int num_chrom);

void density_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
	// variables declaration
	int i,j,k;
	const int nAtoms = residue[atoms[FA->map_par[0].atm].ofres].latm[0] - residue[atoms[FA->map_par[0].atm].ofres].fatm[0] + 1;
	uint maxDensity;
	int mean, stddev;
	int nClusters, nUnclustered;
	float maxDist, minDist;
	float* RMSD;
	double Pi;
	double partition_function;
	ClusterChrom* Chrom;
	ClusterChrom *pChrom, *iChrom, *jChrom;
	Cluster* Clust;
	Cluster* pCluster;

	// dynamic memory allocation
	Chrom = (ClusterChrom*) malloc(num_chrom*sizeof(ClusterChrom));
	RMSD = (float*) malloc(num_chrom*num_chrom*sizeof(float));
	
	//  dynamically allocated memory check-up
	if(!Chrom || !RMSD)
	{
		fprintf(stderr,"ERROR: memory allocation error for clusters\n");
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
			pChrom->isCore = false;
			pChrom->isCenter = false;
			pChrom->Chromosome = &chrom[i];
			pChrom->Cluster = 0;
			pChrom->Density = 0;
			pChrom->CF = 0.0;
			pChrom->DP = NULL;
			pChrom->DPdist = 0.0;
			memset(pChrom->Coord, 0.0, 3*MAX_ATM_HET);
			if(FA->temperature) partition_function += pow( E, ((-1.0) * FA->beta * pChrom->Chromosome->app_evalue) );
		}
	}

	// (0) VERIFICATIONS 
	// Verify that partition_function != NULL
	if(FA->temperature && partition_function == 0.0) 
	{
		fprintf(stderr,"ERROR: The Partition Function is NULL during the clustering step.\n");
		Terminate(2);
	}

	// (1) Build Chromosome Cartesian Coordinates
	for(i = 0; i < num_chrom; ++i)
	{
		pChrom = &Chrom[i];
		if(i+1 == num_chrom) calc_rmsd_chrom(FA,GB,chrom,gen_lim,atoms,residue,cleftgrid,GB->num_genes,i,i, pChrom->Coord, NULL, false);
		else calc_rmsd_chrom(FA,GB,chrom,gen_lim,atoms,residue,cleftgrid,GB->num_genes,i,i+1, pChrom->Coord, (pChrom++)->Coord, false);
	}

	// (2) Build RMSD Matrix
	for(i = 0, Pi = 0.0, iChrom=NULL; i < num_chrom; ++i)
	{
		iChrom = &Chrom[i];
		if(FA->temperature)
		{
			Pi = pow( E, ((-1.0) * FA->beta * iChrom->Chromosome->app_evalue) ) / partition_function;
			iChrom->CF = (double) ( Pi * iChrom->Chromosome->app_evalue) - (FA->temperature * Pi * log(Pi));
		}
		else iChrom->CF = iChrom->Chromosome->app_evalue;
		for(j = i+1, minDist = 0.0; j < num_chrom; ++j)
		{
			jChrom = &Chrom[j];
			for(k = 0; k < nAtoms; ++k) minDist += sqrdist(iChrom->Coord, jChrom->Coord);
			RMSD[K(i,j,num_chrom)] = sqrtf(minDist/(float)num_chrom);
		}
	}

	// (3) Build Local Density Matrix
	for(i = 0; i < num_chrom; ++i)
	{
		iChrom = &Chrom[i];
		// for(j = i+1; j < num_chrom; ++j)
		for(j = 0; j < num_chrom; ++j)
		{
			jChrom	= &Chrom[j];
			if(jChrom != iChrom) iChrom->Density += Chi(RMSD[K(i,j,num_chrom)], FA->cluster_rmsd);
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

	// pChrom (which points to the highest density chrom) should not have pChrom->DP defined
	if(pChrom->DP == NULL)
	{
		pChrom->DPdist = maxDist;
		pChrom->PiDi = pChrom->Density * pChrom->DPdist;
		--k;
	}
	
	// (*) Dealing with multiple Chrom->DP == NULL
	while(k)
	{
		i = 0; pChrom = Chrom; // equivalent to pChrom = &Chrom[0]
		while(pChrom->DP != NULL && pChrom->DPdist > 0.0 && i < num_chrom) {pChrom++; i++;}
		if(pChrom->Density == maxDensity) pChrom->DPdist = maxDist;
		else for(j=0, minDist = FLT_MAX; j < num_chrom; ++j)
		{
			jChrom = &Chrom[j];
			if(jChrom != pChrom && jChrom->Density > pChrom->Density && RMSD[K(i,j,num_chrom)] <= minDist && RMSD[K(i,j,num_chrom)] > 0.0)
			{
				minDist = RMSD[K(i,j,num_chrom)];
				pChrom->DP = jChrom;
				pChrom->DPdist = minDist;
			}
		}
		pChrom->PiDi = pChrom->Density * pChrom->DPdist;
		if(pChrom->PiDi > 0.0) --k;
	}
	

	// (*) Sort Chrom by decreasing PiDi value
	QuickSort_ChromCluster_by_PiDi(Chrom,num_chrom,0,num_chrom-1);
	
	// (6) Identify Cluster Centers
	nUnclustered = num_chrom;
	for(pChrom=NULL, i=0, nClusters=0, stddev=calculate_mean(Chrom, num_chrom), mean=calculate_mean(Chrom, num_chrom); ((&Chrom[nClusters])->PiDi - (&Chrom[nClusters+1])->PiDi) > stddev; ++i)
	{
		pChrom = &Chrom[nClusters];
		
		if(pChrom != NULL)
		{	
			if(pChrom->DP && pChrom->DP->Cluster >= 1) pChrom->Cluster = pChrom->DP->Cluster;
			else
			{
				pChrom->Cluster = (++nClusters);		// cluster assigned
				pChrom->isCenter = true; 				// cluster center assigned
			}
			if(pChrom->Cluster != 0) --nUnclustered;
		}
	}

	// (*) QuickSort by decreasing Density value
	QuickSort_ChromCluster_by_Density(Chrom, num_chrom, 0, num_chrom-1);

	// (7) Clustering Step
	for(i=0, pChrom=NULL; i<num_chrom && nUnclustered > 0; ++i)
	{
		pChrom = &Chrom[i];
		iChrom = &Chrom[i];
		while(pChrom->Cluster <= 0 && pChrom->DP != NULL) pChrom = pChrom->DP;
		if(pChrom->Cluster > 0) iChrom->Cluster = pChrom->Cluster;
		if(iChrom->Cluster > 0 && iChrom!=pChrom) --nUnclustered;
	}

	// (*) At this point, the cluster core vs cluster halo assignation is unuseful if there is a single cluster (unable to separe halore (noise) from core data points in dataset)
	if(nClusters > 1)
	{
		// (8) Assignation of chromosome to its cluster Core/Halo
		// 		1. Find for each Cluster(k): define the border region 
		// 		2. Find for each Cluster(m): the point of highest density(Pb) within the border region
		for(k = 1, pChrom=NULL, maxDensity=0; k <= nClusters; ++k)
		{
			for(i=0, iChrom=Chrom; i<num_chrom; ++i, ++iChrom) if( k == iChrom->Cluster )
			{
				for(j=0, jChrom=Chrom; j<num_chrom; ++j, ++jChrom) if(jChrom->Cluster > 0 && jChrom->Cluster != k && RMSD[K(iChrom->index, jChrom->index, num_chrom)] <= FA->cluster_rmsd && iChrom->Density > maxDensity)
				{
					pChrom = iChrom;
					maxDensity = pChrom->Density;
				}
			}
			for(i=0, iChrom=Chrom; i<num_chrom;++i, ++iChrom) if( k == iChrom->Cluster && pChrom != NULL)
			{
				if(iChrom->Density >= pChrom->Density) pChrom->isCore = true;
			}
		}

		// (9) Cluster Creation
		Clust = (Cluster*) malloc(nClusters*sizeof(Cluster));
		//  dynamically allocated memory check-up
		if(Clust == NULL) 
		{
			fprintf(stderr,"ERROR: memory allocation error for clusters\n");
			Terminate(2);
		}
		for(k = 1, pCluster=Cluster; k <= nClusters && pCluster != NULL; ++k, ++pCluster)
		{
			if(pCluster != &Clust[k]) break; 	 // pCluster should always equal &Clust[k]
			// initializing the Cluster element
			pCluster->ID = k;
			pCluster->Frequency = 0;
			pCluster->totCF = 0.0;
			pCluster->BestCF = NULL;
			pCluster->Center = NULL;

			for(pChrom=Chrom, i=0; i<num_chrom; ++i, ++pChrom) if(pChrom->isCore && pChrom->Cluster == k)
			{
				pCluster->totCF += pChrom->CF;
				++(pCluster->Frequency);
				if(pChrom->isCenter) pCluster->Center = pChrom;
				if( (pCluster->BestCF == NULL) || (pChrom->Chromosome->app_evalue < pCluster->BestCF->Chromosome->app_evalue) ) pCluster->BestCF = pChrom;
			}
		}
	}
	// else if(nClusters == 1)
	// {

	// }

	// (*) Printing informations
    for(i=0,j=0;i<num_chrom;++i)
    {
        printf("Add:%p\tDensity:%d\tDistance:%g\tCluster:%d\tDP:%p\tPiDi:%g\tCF:%g\tisCore:%s\n",&Chrom[i],(&Chrom[i])->Density,(&Chrom[i])->DPdist,(&Chrom[i])->Cluster,(&Chrom[i])->DP,(&Chrom[i])->PiDi,(&Chrom[i])->CF,(&Chrom[i])->isCore ? "true" : "false");
        if((&Chrom[i])->DP != NULL && (&Chrom[i])->Density > (&Chrom[i])->DP->Density) j++;
    }
    printf("There is %d chromosomes which DP are of lower Density.\nMean:%g\tStdDev:%g\n",j,calculate_mean(Chrom, num_chrom),calculate_stddev(Chrom, num_chrom));
	
	// (*) Memory deallocation
	if(Chrom) { free(Chrom); Chrom=NULL; } 
	if(RMSD)  { free(RMSD); RMSD=NULL; }

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
	ClusterChrom* ChromT = NULL;
	
	ChromT = (ClusterChrom*)malloc(sizeof(ClusterChrom));
	if(ChromT == NULL)
	{
		fprintf(stderr,"ERROR: memory allocation error for clusters\n");
		Terminate(2);
	}

	ChromT->index = ChromX->index; ChromX->index = ChromY->index; ChromY->index = ChromT->index;
	ChromT->isCore = ChromX->isCore; ChromX->isCore = ChromY->isCore; ChromY->isCore = ChromT->isCore;
	ChromT->isCenter = ChromX->isCenter; ChromX->isCenter = ChromY->isCenter; ChromY->isCenter = ChromT->isCenter;
	ChromT->Chromosome = ChromX->Chromosome; ChromX->Chromosome = ChromY->Chromosome; ChromY->Chromosome = ChromT->Chromosome;
	ChromT->Cluster = ChromX->Cluster; ChromX->Cluster = ChromY->Cluster; ChromY->Cluster = ChromT->Cluster;
	ChromT->CF = ChromX->CF; ChromX->CF = ChromY->CF; ChromY->CF = ChromT->CF;
	ChromT->Density = ChromX->Density; ChromX->Density = ChromY->Density; ChromY->Density = ChromT->Density;
	ChromT->PiDi = ChromX->PiDi; ChromX->PiDi = ChromY->PiDi; ChromY->PiDi = ChromT->PiDi;
	ChromT->DPdist = ChromX->DPdist; ChromX->DPdist = ChromY->DPdist; ChromY->DPdist = ChromT->DPdist;
	
	for(i=0, pChrom=NULL;i<num_chrom;++i)
	{
		pChrom = &Chrom[i];
		// Check these conditions
		if(/*pChrom != ChromX && pChrom != ChromY && */pChrom->DP == ChromX) pChrom->DP = ChromY;
		else if(/*pChrom != ChromX && pChrom != ChromY && */pChrom->DP == ChromY) pChrom->DP = ChromX;
	}
	
	ChromT->DP = ChromX->DP; ChromX->DP = ChromY->DP; ChromY->DP = ChromT->DP;
	
	for(i=0;i<3*MAX_ATM_HET;++i)
	{
		ChromT->Coord[i] = ChromX->Coord[i]; ChromX->Coord[i] = ChromY->Coord[i]; ChromY->Coord[i] = ChromT->Coord[i];
	}
	
	if(ChromT) { free(ChromT); ChromT=NULL; }
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