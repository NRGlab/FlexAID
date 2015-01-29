#include "gaboom.h"
#include "boinc.h"
#include "float.h"

// the following macro is used to retrieve the appropriate index (k) or indices (i,j) of a 1D array instead of an equivalent 2D matrix and vice-versa
#define K(i, j, n) 	( (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1 )
#define I(k,n) 		( n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5))
#define J(k,n) 		( k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2 )
// Chi(x) function as defined in :
//   Science 344(6191):1492-1496, 2014, Eq. (1)
#define Chi(a,d) 	( ((a-d) < 0.0) ? 1 : 0 )

// Function prototypes
void calculate_stddev(	float* PiDi, int num_chrom);
void QuickSort_Density(	chromosome* chrom, double* appCF, float* di, int* density, float* PiDi, int* assigned_cluster, int beg, int end);
void QuickSort_PiDi(	chromosome* chrom, double* appCF, float* di, int* density, float* PiDi, int* assigned_cluster, int beg, int end);
void swap_elements(		chromosome* CHROMx, double* appCFx, float* DISTx, int* DENSITYx, float* PiDix, int* CLUSx,
						chromosome* CHROMy, double* appCFy, float* DISTy, int* DENSITYy, float* PiDiy, int* CLUSy);
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
	int dim = ( ( num_chrom * (num_chrom + 1) ) / 2 - num_chrom );	/* 1D array dimension needed to emulate storage of a symmetric 2D matrix */
	int i = 0, j = 0, k = 0;										/* indexing iterators */
	int nClusters = 0;												/* number of defined clusters */
	int nOutliers = 0;												/* number of defined outliers (chrom that has not been assigned to a cluster) */
	int nUnclustered = 0;
	float stddev = 0.0f;
	float* rmsd_matrix = NULL;										/* Step (1) */
	int* density_matrix = NULL;										/* Step (2) */ 
	float* distance_matrix = NULL;									/* Step (3) */ 
	int* assigned_cluster = NULL;									/* Assigned Cluster (0:unclustered, -1:outlier) */
	int* nearest_center = NULL;
	double* appCF = NULL;
	float* PiDi = NULL;
	// memory allocation
	rmsd_matrix = (float*)malloc(dim*sizeof(float)); 			// dim * sizeof is used to simulate 2D array of chrom*chrom (real size of dim = chrom(chrom+1)/2)
	density_matrix = (int*)malloc(num_chrom*sizeof(int));		// num_chrom * sizeof is used to build a 1D array
	distance_matrix = (float*)malloc(num_chrom*sizeof(float));	// num_chrom * sizeof is used to build a 1D array
	assigned_cluster = (int*)malloc(num_chrom*sizeof(int));		/* Assigned Cluster (0:unclustered, -1:outlier) is 1D array sizeof(num_chrom*(int)) */
	nearest_center = (int*)malloc(num_chrom*sizeof(int));		/* Assigned Cluster (0:unclustered, -1:outlier) is 1D array sizeof(num_chrom*(int)) */
	appCF = (double*)malloc(num_chrom*sizeof(double));			/* appCF */
	PiDi = (float*)malloc(num_chrom*sizeof(float));
	// memory check
	if(!rmsd_matrix || !density_matrix || !distance_matrix || !appCF || !assigned_cluster || !PiDi || nearest_center)
	{
		fprintf(stderr,"ERROR: memory allocation error for clusters\n");
		Terminate(2);
	}

	// matrices initialization
	for(i = 0; i <= dim; ++i) rmsd_matrix[i] = 0.0f;
	for(i = 0; i <= num_chrom; ++i)
	{
		density_matrix[i] = 0;
		distance_matrix[i] = 0.0f;
		assigned_cluster[i] = 0;
		appCF[i] = 0.0;
		// line below uses of this loop to iterate through all chrom to calculate the partition_function of the population
		if(FA->temperature) partition_function += pow( E, ((-1.0) * FA->beta * chrom[i].app_evalue) );
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
		for(j = 0; j < num_chrom; ++j) 
			if(i != j) rmsd_matrix[ K(i,j,num_chrom) ] = calc_rmsd_chrom(FA,GB,chrom,gen_lim,atoms,residue,cleftgrid,GB->num_genes,i,j);
	}

	
	// (2) build local density matrix
	for(i = 0; i < num_chrom; ++i) 
		for(j = 0; j < num_chrom; ++j) 
			if(i != j) density_matrix[i] += Chi( rmsd_matrix[K(i,j,num_chrom)], FA->cluster_rmsd );

	
	// (3) build distance_matrix from rmsd_matrix && density_matrix
	for(i = 0; i < num_chrom; ++i)
	{
		for(j = 0, distance = FLT_MAX; j < num_chrom; ++j) if( i != j && density_matrix[j] > density_matrix[i] && rmsd_matrix[K(i,j,num_chrom)] < distance) distance = rmsd_matrix[K(i,j,num_chrom)];
		distance_matrix[i] = distance;
		assigned_cluster[i] = j;
	}
	if(rmsd_matrix != NULL) free(rmsd_matrix); // Freeing rmsd_matrix before sorting clusters, as pairwise RMSD is no longer needed.

	// Setting a ∂i value to the chrom with highest density_matrix value (distance will have the max)
	for(i=0, j=0, k=0; i<num_chrom; ++i) if(distance_matrix[i] > distance_matrix[j]) 
	{
		if(distance_matrix[i] < FLT_MAX	)
		{
			j = i;
			PiDi[i] = distance_matrix[i] * density_matrix[i]; // Profit of this loop to compute PiDi !
		}
		else k = i;
	}
	distance_matrix[k] = distance_matrix[j];
	PiDi[k] = distance_matrix[k] * density_matrix[k]; // Filling out 

	// (4) Order chromosomes in decreasing ∂i values (which are stored in density_matrix[])
	QuickSort_PiDi(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, 0, num_chrom-1);

	
	// (5) Count Clusters and assign cluster centers
	stddev = calculate_stddev(PiDi, num_chrom); nClusters = 0;
	while( (PiDi[nClusters] - PiDi[nClusters+1]) > 2*stddev) ++nClusters; // Check the stop condition on this loop

	
	// (6) Identify Cluster Centers

	
	// freeing dynamically allocated memory
	if(appCF != NULL) free(appCF);
	if(assigned_cluster != NULL) free(assigned_cluster);
	if(density_matrix != NULL) free(density_matrix);
	if(distance_matrix != NULL) free(distance_matrix);
	if(PiDi != NULL) free(PiDi);
	if(nearest_center) free(nearest_center);
}

void QuickSort_PiDi(chromosome* chrom, double* appCF, float* distance_matrix, int* density_matrix, float* PiDi, int* assigned_cluster, int beg, int end)
{
	int l,r,p;
	int value;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		value = PiDi[p];

		while(1)
		{
			while( (l<=r) && QS_DSC(density_matrix[l],value) <= 0 ) ++l;
			while( (l<=r) && QS_DSC(density_matrix[r],value)  > 0 ) --r;
			
			if (l > r) break;

			swap_elements(	&chrom[p],&appCF[p],&distance_matrix[p],&density_matrix[p],&PiDi[p], &assigned_cluster[p],
							&chrom[r],&appCF[r],&distance_matrix[r],&density_matrix[r],&PiDi[r], &assigned_cluster[r] );
			--r;

			if( (r-beg) < (end-l) )
			{
				QuickSort_PiDi(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, beg, r);
				beg = l;
			}
			else
			{
				QuickSort_PiDi(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, l, end);
				end = r;
			}
		}
	}
}
// this function is used to reorder cluster by decreasing density
void QuickSort_Density(chromosome* chrom, double* appCF, float* distance_matrix, int* density_matrix, float* PiDi, int* assigned_cluster, int beg, int end)
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

			swap_elements(	&chrom[p],&appCF[p],&distance_matrix[p],&density_matrix[p],&PiDi[p], &assigned_cluster[p],
							&chrom[r],&appCF[r],&distance_matrix[r],&density_matrix[r],&PiDi[r], &assigned_cluster[r] );
			--r;

			if( (r-beg) < (end-l) )
			{
				QuickSort_Density(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, beg, r);
				beg = l;
			}
			else
			{
				QuickSort_Density(chrom, appCF, distance_matrix, density_matrix, PiDi, assigned_cluster, l, end);
				end = r;
			}
		}
	}
}
void swap_elements(	chromosome* CHROMx, double* appCFx, float* DISTx, int* DENSITYx, float* PiDix,int* CLUSx,
					chromosome* CHROMy, double* appCFy, float* DISTy, int* DENSITYy, float* PiDiy,int* CLUSy )
{
	// temporary variables 
	chromosome CHROMt;
	double appCFt = 0.0;
	float RMSDt = 0.0f, DISTt = 0.0f, PiDit = 0.0f;
	int DENSITYt = 0, CLUSt = 0;
	CHROMt = *CHROMx; *CHROMx = *CHROMy; *CHROMy = CHROMt;
	DISTt = *DISTx; *DISTx = *DISTy; *DISTy = DISTt;
	PiDit = *PiDix; *PiDix = *PiDiy; *PiDiy = PiDit;
	DENSITYt = *DENSITYx; *DENSITYx = *DENSITYy; *DENSITYy = DENSITYt;
	CLUSt = *CLUSx; *CLUSx = *CLUSy; *CLUSy = CLUSt;
}

float calculate_stddev(float* PiDi, int num_chrom)
{
	int i;
	float mean = 0.0f;
	float sqrtot = 0.0f;
	for(i = 0; i < num_chrom; ++i) mean += PiDi[i];
	mean /= num_chrom;
	for(i = 0; i < num_chrom; ++i) sqrtot += pow((mean - PiDi[i]),2);
	return sqrtf(sqrtot);
}
