#include "gaboom.h"
#include "boinc.h"

// Chi(x) function as defined in :
//   Science 344(6191):1492-1496, 2014, Eq. (1)
#define Chi(a,d) 	( ((a-d) < 0.0) ? 1 : 0 )
#define K(i,j,n) ( (i < j) ? (i*n+j) : (j*n+i) )
// Function prototypes
float calculate_stddev(float* PiDi, int num_chrom);
void QuickSort_Density(chromosome* chrom, double* appCF, float* di, int* density, float* PiDi, int* assigned_cluster, int* nearest_center, int beg, int end);
void QuickSort_PiDi(chromosome* chrom, double* appCF, float* di, int* density, float* PiDi, int* assigned_cluster, int* nearest_center, int beg, int end);
void swap_elements(	chromosome* CHROMx, double* appCFx, float* DISTx, int* DENSITYx, float* PiDix, int* CLUSx, int* CENx,
						chromosome* CHROMy, double* appCFy, float* DISTy, int* DENSITYy, float* PiDiy, int* CLUSy, int* CENy);
// Clustering structures
struct cluster{
	int id;
	int nChrom;
	double totCF;
	float refRMSD;
	chromosome* chromosomes;
	chromosome* representative;
	chromosome* best;
};
typedef struct cluster Cluster;

struct clusters{
	int nClusters;
	int nOutliers;

	Cluster* first,last;
	Cluster* next, previous;
	chromosome* outliers; 
};
typedef struct clusters Clusters;


// Main Function
/*
	density_cluster(**args) is ...  
*/
void density_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
	// The following steps will be executed prior to assign clusters:
	/*	
		(1) Build RSMD matrix between each chrom [dij]
		(2) Build Local Density matrix of each chrom (Eq. 1) [pij] : pi = ∑jCHI(dij - d_cutoff)
		(3) Build Distance from Density Centers matrix (Eq. 2) [∂i] = min(dij) for pj > pi
		(4) Order chrom in decreasing ∂i : ∂i is much larger for chrom which are local/global maxima in density
		(5) Count Clusters (count centers) Cluster centers are recognized as chrom with anomalously large ∂i values
		(6) 
	*/

	// variables initialization
	bool Hungarian = false; 										/* Boolean indicating whether to use symmetry corrected RMSD or not */
	// cfstr cf;														/* Complementarity Function value */
	double Pi = 0.0;												/* variable that will be iteratively calculated */ 
	double partition_function = 0.0;								/* will represent the partition_function of chromosomes CF values */ 
	float distance = FLT_MAX;										/* float variable that will be used as temporary memory for RMSD values (float) */
	int i = 0, j = 0, k = 0;										/* indexing iterators */
	int nClusters = 0;												/* number of defined clusters */
	int nOutliers = 0;												/* number of defined outliers (chrom that has not been assigned to a cluster) */
	int nUnclustered = 0;
	float stddev = 0.0f;
	int dim = num_chrom*num_chrom;
	float* rmsd_matrix = NULL;										/* Step (1) */
	int* density_matrix = NULL;										/* Step (2) */ 
	float* distance_matrix = NULL;									/* Step (3) */ 
	int* assigned_cluster = NULL;									/* Assigned Cluster (0:unclustered, -1:outlier) */
	int* nearest_center = NULL;
	double* appCF = NULL;
	float* PiDi = NULL;
	
	// memory allocation
	density_matrix = (int*)malloc(num_chrom*sizeof(int));		// num_chrom * sizeof is used to build a 1D array
	distance_matrix = (float*)malloc(num_chrom*sizeof(float));	// num_chrom * sizeof is used to build a 1D array
	assigned_cluster = (int*)malloc(num_chrom*sizeof(int));		/* Assigned Cluster (0:unclustered, -1:outlier) is 1D array sizeof(num_chrom*(int)) */
	nearest_center = (int*)malloc(num_chrom*sizeof(int));		/* Assigned Cluster (0:unclustered, -1:outlier) is 1D array sizeof(num_chrom*(int)) */
	appCF = (double*)malloc(num_chrom*sizeof(double));			/* appCF */
	PiDi = (float*)malloc(num_chrom*sizeof(float));
	rmsd_matrix = (float*)malloc(dim*sizeof(float));
	
	// memory check
	if(!rmsd_matrix || !density_matrix || !distance_matrix || !appCF || !assigned_cluster || !PiDi || !nearest_center)
	{
		fprintf(stderr,"ERROR: memory allocation error for clusters\n");
		Terminate(2);
	}
	
	// matrices initialization
	memset(rmsd_matrix, 0.0, dim);
	for(i = 0; i < num_chrom; ++i)
	{
		density_matrix[i] = 0;
		distance_matrix[i] = 0.0;
		assigned_cluster[i] = 0;
		nearest_center[i] = 0;
		appCF[i] = 0.0;
		// line below uses this loop to iterate through all chrom to calculate the partition_function of the population
		if(FA->temperature) partition_function += pow( E, ((-1.0) * FA->beta * chrom[i].app_evalue) );
	}

	// (0) VERIFICATIONS 
	// Verify that partition_function != NULL
	if(FA->temperature && partition_function == 0.0) 
	{
		fprintf(stderr,"ERROR: The Partition Function is NULL in the clustering step.\n");
		Terminate(2);
	}

	// (1) build RMSD matrix between each chrom
	for(i = 0; i < num_chrom; ++i)
	{
		// if FA->temperature: the conformational entropy will be calculated in the partition function
		if(FA->temperature)
		{
			Pi = pow( E, ((-1.0) * FA->beta * chrom[i].app_evalue) ) / partition_function;
			appCF[i] = (double)( ( Pi * chrom[i].app_evalue) - (FA->temperature * Pi * log(Pi)) );
		}
		// else: the normal apparent CF value will be used
		else appCF[i] = chrom[i].app_evalue;
		// INNER LOOP: calculate pairwise RMSD values between each chrom and stores if in rmsd_matrix[]
		for(j = i+1; j < num_chrom; ++j) rmsd_matrix[ K(i,j,num_chrom) ] = calc_rmsd_chrom(FA,GB,chrom,gen_lim,atoms,residue,cleftgrid,GB->num_genes,i,j);
	}

	
	// (2) build local density matrix
	for(i = 0; i < num_chrom; ++i)
    {
        for(j=i+1; j < num_chrom; ++j) density_matrix[i] += Chi( rmsd_matrix[K(i,j,num_chrom)], FA->cluster_rmsd );
    }
	
	// (3) build distance_matrix from rmsd_matrix && density_matrix
    for(i = 0; i < num_chrom; ++i)
    {
        for(j = 0, distance = FLT_MAX; j < num_chrom; ++j)
        {
            if(i!=j && density_matrix[j] > density_matrix[i] && rmsd_matrix[K(i,j,num_chrom)] < distance )
            {
                distance = rmsd_matrix[K(i,j,num_chrom)];
                nearest_center[i] = j;
                distance_matrix[i] = distance;
            }
        }
    }
	// Setting a ∂i value to the chrom with highest density_matrix value (distance will have the maximal distance)
	// Profit of this loop to compute Gamma = PiDi !
	for(i=0, j=0; i<num_chrom; ++i)
	{
		(distance_matrix[i] > distance_matrix[j]) ? (j = i) : (PiDi[i] = distance_matrix[i] * density_matrix[i]);
	}

    // (4) Order chromosomes in decreasing Pi∂i values (which are stored in density_matrix[])
	// if(rmsd_matrix != NULL) free(rmsd_matrix); // Freeing rmsd_matrix before sorting clusters, as pairwise RMSD is no longer needed.
    // QuickSort_Density(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, nearest_center, 0, num_chrom-1);
    // (*) Printing
    for(i = 0; i<num_chrom;++i) printf("i:%3d\tdensity:%3d\tdistance:%g\tCF:%g\tnearest center:%d\tPiDi:%g\n",i,density_matrix[i], distance_matrix[i], appCF[i], nearest_center[i], PiDi[i]);

    // (5) Count Clusters and Identify Cluster Centers
	stddev = calculate_stddev(PiDi, num_chrom); nClusters = 1; nUnclustered = num_chrom;
    while( (PiDi[nClusters] - PiDi[nClusters+1]) > 2*stddev )  // while the difference in Pi*Di is anomalously large (please Science define the fuck out of this)
	{
		assigned_cluster[nClusters] = nClusters;
		++nClusters;
        --nUnclustered;
	}	
	
	// (6) Clustering step!
	QuickSort_Density(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, nearest_center, 0, num_chrom-1);
	while(nUnclustered)
	{
		for(i = 0; i < num_chrom; ++i) if(assigned_cluster[i] < 1) break; // identify the next unclustered chrom
		
        // Cluster HERE
	}

	// freeing dynamically allocated memory
	if(appCF != NULL) 			free(appCF);
	if(assigned_cluster != NULL)free(assigned_cluster);
	if(density_matrix != NULL) 	free(density_matrix);
	if(distance_matrix != NULL) free(distance_matrix);
	if(PiDi != NULL) 			free(PiDi);
	if(nearest_center != NULL)	free(nearest_center);
}

void QuickSort_PiDi(chromosome* chrom, double* appCF, float* distance_matrix, int* density_matrix, float* PiDi, int* assigned_cluster, int *nearest_center, int beg, int end)
{
	int l,r,p;
	float value;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		value = PiDi[p];

		while(1)
		{
			while( (l<=r) && QS_DSC(PiDi[l],value) <= 0 ) ++l;
			while( (l<=r) && QS_DSC(PiDi[r],value)  > 0 ) --r;
			
			if (l > r) break;

			swap_elements(	&chrom[l], &appCF[l], &distance_matrix[l], &density_matrix[l], &PiDi[l], &assigned_cluster[l], &nearest_center[l],
							&chrom[r], &appCF[r], &distance_matrix[r], &density_matrix[r], &PiDi[r], &assigned_cluster[r], &nearest_center[r]);
			if (p == r) p=l;
			++l;--r;
		}
		swap_elements(	&chrom[p], &appCF[p], &distance_matrix[p], &density_matrix[p], &PiDi[p], &assigned_cluster[p], &nearest_center[p],
						&chrom[r], &appCF[r], &distance_matrix[r], &density_matrix[r], &PiDi[r], &assigned_cluster[r], &nearest_center[r]);
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_PiDi(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, nearest_center, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_PiDi(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, nearest_center, l, end);
			end = r;
		}
	}
}
// this function is used to reorder cluster by decreasing density
void QuickSort_Density(chromosome* chrom, double* appCF, float* distance_matrix, int* density_matrix, float* PiDi, int* assigned_cluster, int* nearest_center, int beg, int end)
{
	int l,r,p;
	int value;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		value = density_matrix[p];

		while(1)
		{
			while( (l<=r) && QS_DSC(density_matrix[l],value) <= 0 ) ++l;
			while( (l<=r) && QS_DSC(density_matrix[r],value)  > 0 ) --r;
			
			if (l > r) break;

			swap_elements(	&chrom[l],&appCF[l],&distance_matrix[l],&density_matrix[l],&PiDi[l], &assigned_cluster[l], &nearest_center[l],
							&chrom[r],&appCF[r],&distance_matrix[r],&density_matrix[r],&PiDi[r], &assigned_cluster[r], &nearest_center[r]	);
			if (p == r) p=l;
			++l;--r;
		}
		swap_elements(	&chrom[p],&appCF[p],&distance_matrix[p],&density_matrix[p],&PiDi[p], &assigned_cluster[p], &nearest_center[p],
						&chrom[r],&appCF[r],&distance_matrix[r],&density_matrix[r],&PiDi[r], &assigned_cluster[r], &nearest_center[r]	);
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_Density(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, nearest_center, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_Density(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, nearest_center, l, end);
			end = r;
		}
	}
}
void swap_elements(	chromosome* CHROMx, double* appCFx, float* DISTx, int* DENSITYx, float* PiDix, int* CLUSx, int*CENx,
					chromosome* CHROMy, double* appCFy, float* DISTy, int* DENSITYy, float* PiDiy, int* CLUSy, int*CENy)
{
	// temporary variables 
	chromosome CHROMt;
	double appCFt = 0.0;
	float DISTt = 0.0f, PiDit = 0.0f;
	int DENSITYt = 0, CLUSt = 0, CENt = 0;
	appCFt = *appCFx; *appCFx = *appCFy; *appCFy = appCFt;
	CHROMt = *CHROMx; *CHROMx = *CHROMy; *CHROMy = CHROMt;
	DISTt = *DISTx; *DISTx = *DISTy; *DISTy = DISTt;
	PiDit = *PiDix; *PiDix = *PiDiy; *PiDiy = PiDit;
	DENSITYt = *DENSITYx; *DENSITYx = *DENSITYy; *DENSITYy = DENSITYt;
	CLUSt = *CLUSx; *CLUSx = *CLUSy; *CLUSy = CLUSt;
	CENt = *CENx; *CENx = *CENy; *CENy = CENt;
}

float calculate_stddev(float* PiDi, int num_chrom)
{
	int i;
	float mean = 0.0f;
	float sqrtot = 0.0f;
    // printf("float max = %g\n", FLT_MAX);
    for(i = 0; i < num_chrom; ++i)
    {
        mean += PiDi[i];
    }
	mean /= num_chrom;
    // printf("mean PiDi = %g\n", mean);
	for(i = 0; i < num_chrom; ++i)
    {
        sqrtot += pow( (mean - PiDi[i]),2.0 );
    }
    printf("mean:%g\tstddev:%g\n",mean,sqrtf(sqrtot/(float)num_chrom));
	return sqrtf(sqrtot/(float)num_chrom);
}
