#include "gaboom.h"
#include "boinc.h"

// Chi(x) function as defined in :
//   Science 344(6191):1492-1496, 2014, Eq. (1)
#define Chi(a,d) 	( ((a-d) < 0.0) ? 1 : 0 )
#define K(i,j,n) ( (i < j) ? (i*n+j) : (j*n+i) )
// Function prototypes
float calculate_stddev(float* PiDi, int num_chrom);
float calculate_mean(float* PiDi, int num_chrom);
void  QuickSort_Density(chromosome* chrom, double* appCF, float* di, int* density, int* assigned_cluster, int* nearest_center, float* PiDi, float* coord, float* rmsd, int num_chrom, int beg, int end);
void  QuickSort_PiDi(chromosome* chrom, double* appCF, float* di, int* density, int* assigned_cluster, int* nearest_center, float* PiDi, float* coord, float* rmsd, int num_chrom, int beg, int end);
void  swap_elements(chromosome* CHROMx, double* appCFx, float* DISTx, int* DENSITYx, int* CLUSx, int* CENx, float* PiDix, float* COORDx, float* RMSDx, chromosome* CHROMy, double* appCFy, float* DISTy, int* DENSITYy, int* CLUSy, int* CENy, float* PiDiy, float* COORDy, float* RMSDy, int num_chrom);
void  assign_cluster_from_density_neighborhood(int i, int* nUnclustered, int* nOutliers, int* density_matrix, float* distance_matrix, int* nearest_center, int* assigned_cluster);

// Clustering structures
struct Cluster{
	int id;
	int size;
	double totCF;
	float refRMSD;
	chromosome* chromosomes;
	chromosome* representative;
	chromosome* best;
};
typedef struct Cluster Cluster;

struct Clusters{
	int nClusters;
	int nOutliers;

	Cluster* first,last;
	Cluster* next, previous;
	chromosome* outliers; 
};
typedef struct Clusters Clusters;


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
	int i, j, k, l,m;													/* indexing iterators */
	int max_density = -1;
	int nAtoms;
	int nClusters = 0;												/* number of defined clusters */
	int nOutliers = 0;												/* number of defined outliers (chrom that has not been assigned to a cluster) */
	int nUnclustered = 0;
	float stddev = 0.0f, mean = 0.0f;
    int dim = num_chrom*num_chrom;
	float* rmsd_matrix = NULL;										/* Step (1) */
	int* density_matrix = NULL;										/* Step (2) */ 
	float* distance_matrix = NULL;									/* Step (3) */ 
	int* assigned_cluster = NULL;									/* Assigned Cluster (0:unclustered, -1:outlier) */
	int* nearest_center = NULL;
	double* appCF = NULL;
	float* PiDi = NULL;
	float* cartesian_coord = NULL;
	// memory allocation
	cartesian_coord = (float*)malloc(num_chrom*3*MAX_ATM_HET*sizeof(float));
	density_matrix = (int*)malloc(num_chrom*sizeof(int));		// num_chrom * sizeof is used to build a 1D array
	distance_matrix = (float*)malloc(num_chrom*sizeof(float));	// num_chrom * sizeof is used to build a 1D array
	assigned_cluster = (int*)malloc(num_chrom*sizeof(int));		/* Assigned Cluster (0:unclustered, -1:outlier) is 1D array sizeof(num_chrom*(int)) */
	nearest_center = (int*)malloc(num_chrom*sizeof(int));		/* Assigned Cluster (0:unclustered, -1:outlier) is 1D array sizeof(num_chrom*(int)) */
	appCF = (double*)malloc(num_chrom*sizeof(double));			/* appCF */
	PiDi = (float*)malloc(num_chrom*sizeof(float));
	rmsd_matrix = (float*)malloc(dim*sizeof(float));
	
	// memory check
	if(!rmsd_matrix || !density_matrix || !distance_matrix || !appCF || !assigned_cluster || !PiDi || !nearest_center || !cartesian_coord)
	{
		fprintf(stderr,"ERROR: memory allocation error for clusters\n");
		Terminate(2);
	}
	l=atoms[FA->map_par[0].atm].ofres;
	nAtoms=residue[l].latm[0] - residue[l].fatm[0] + 1;

	// matrices initialization
	memset(rmsd_matrix, 0.0, dim);
	for(i = 0; i < num_chrom; ++i)
	{
		density_matrix[i] = 0;
		distance_matrix[i] = 0.0;
		assigned_cluster[i] = 0;
		nearest_center[i] = -1;
		appCF[i] = 0.0;
		// line below uses this loop to iterate through all chrom to calculate the partition_function of the population
		if(FA->temperature) partition_function += pow( E, ((-1.0) * FA->beta * chrom[i].app_evalue) );
	}

	// (0) VERIFICATIONS 
	// Verify that partition_function != NULL
	if(FA->temperature && partition_function == 0.0) 
	{
		fprintf(stderr,"ERROR: The Partition Function is NULL during the clustering step.\n");
		Terminate(2);
	}
	for(i=0;i<num_chrom;i+=2)
	{
		if(i+1 == num_chrom) calc_rmsd_chrom(FA,GB,chrom,gen_lim,atoms,residue,cleftgrid,GB->num_genes,i,i, &cartesian_coord[3*i*MAX_ATM_HET], NULL, false);
		else calc_rmsd_chrom(FA,GB,chrom,gen_lim,atoms,residue,cleftgrid,GB->num_genes,i,i+1, &cartesian_coord[3*i*MAX_ATM_HET], &cartesian_coord[3*(i+1)*MAX_ATM_HET], false);
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
		for(j = i+1; j < num_chrom; ++j)
		{
			for(k = 0, distance = 0.0; k < nAtoms; ++k)
            {
                distance += sqrdist(&cartesian_coord[i*3*MAX_ATM_HET+k],&cartesian_coord[j*3*MAX_ATM_HET+k]);
            }
			rmsd_matrix[ K(i,j,num_chrom) ] = sqrtf(distance/(float)nAtoms);
		}
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
            if(i!=j && density_matrix[j] > density_matrix[i] && rmsd_matrix[K(i,j,num_chrom)] < distance && rmsd_matrix[K(i,j,num_chrom)] > 0.0)
            {
                distance = rmsd_matrix[K(i,j,num_chrom)];
                nearest_center[i] = j;
                distance_matrix[i] = distance;
            }
        }
       for(j = num_chrom-1, distance = FLT_MAX; j > 0; --j)
       {
           if(i!=j && density_matrix[j] > density_matrix[i] && rmsd_matrix[K(i,j,num_chrom)] < distance && rmsd_matrix[K(i,j,num_chrom)] > 0.0)
           {
               distance = rmsd_matrix[K(i,j,num_chrom)];
               nearest_center[i] = j;
               distance_matrix[i] = distance;
           }
       }
    }
	
	// Getting the maximal Di value
	// Setting a Di value if nearest_center != -1 
	// Profit of this loop to compute PiDi[i] !
	for(i=0, j=0, k=0, m=0, max_density = -1; i<num_chrom; ++i)
	{
		// save highest_distance found
		if(density_matrix[i] > density_matrix[m])
		{
			m = i;
			max_density = density_matrix[m];
		}
		if(i!=j && distance_matrix[i] > distance_matrix[j] ) j = i;
		// In case distance haven't yet been defined
		if(i!=j && !distance_matrix[i] && nearest_center[i] > -1 && density_matrix[nearest_center[i]] > density_matrix[i]) distance_matrix[i] = rmsd_matrix[K(i,nearest_center[i],num_chrom)];		
		// k will serve to count the number of !distance_matrix[i] && nearest_center[i] == -1 occurences
		if(nearest_center[i] == -1) ++k;
		// Compute PiDi[i]
		PiDi[i] = distance_matrix[i] * density_matrix[i];
	}
	
	// retriving the maximal density point (k is mostly supposed to equal 1 but can be >1 when multiple points are of equal highest density values)
	if(nearest_center[m] == -1) {distance_matrix[m] = distance_matrix[j]; --k; PiDi[m] = distance_matrix[m] * density_matrix[m];}
	while(k)
	{	
		i = 0;
		while(nearest_center[i] != -1) ++i;			// This loop increments i until we find a point which has no nearest higher density value points (nearest_center[i] == -1)
		if(density_matrix[i] == max_density) distance_matrix[i] = distance_matrix[j]; 			// j memorizes the position of highest value withn distance_matrix
        else for(l=0, distance=FLT_MAX;l<num_chrom;++l) if(l != i && density_matrix[l] > density_matrix[i] && rmsd_matrix[K(i,l,num_chrom)] <= distance && rmsd_matrix[K(i,l,num_chrom)] > 0.0)
        {
        	distance =rmsd_matrix[K(i,l,num_chrom)];
            nearest_center[i] = l;
            distance_matrix[i] = distance;  
        }
        PiDi[i] = distance_matrix[i] * density_matrix[i];   // refreses the PiDi value after modifying distance_matrix[i] value
		--k;
	}
    
    // (4) Order chromosomes in decreasing Pi∂i values (which are stored in PiDi[])
    QuickSort_PiDi(chrom, appCF, distance_matrix, density_matrix, assigned_cluster, nearest_center, PiDi, cartesian_coord, rmsd_matrix, num_chrom, 0, num_chrom-1);
    
    // (5) Count Clusters
    nUnclustered = num_chrom; // total number of unclustered chrom as of now
    for(nClusters=0,stddev=calculate_stddev(PiDi, num_chrom), mean=calculate_mean(PiDi, num_chrom); PiDi[nClusters]-PiDi[nClusters+1] > mean+stddev; ++nClusters)
    {
    	assigned_cluster[nClusters] = nClusters+1;
    	--nUnclustered;
    }
    
    // (6) Identify cluster centers

    // (7) Clustering Step

    // (*) Density sorting
	// QuickSort_Density(chrom, appCF, distance_matrix, density_matrix, assigned_cluster, nearest_center, PiDi, cartesian_coord, rmsd_matrix, num_chrom, 0, num_chrom-1);
    // (*) Printing
    printf("number of clusters:%d\n",nClusters);
    for(i = 0; i<num_chrom;++i) printf("i:%3d\tdensity:%3d\tdistance:%1.4g\tCF:%4.4g\tnearest center:%d\tPiDi:%g\tCluster:%d\n",i,density_matrix[i], distance_matrix[i], appCF[i], nearest_center[i], PiDi[i], assigned_cluster[i]);

	// freeing dynamically allocated memory
	if(appCF != NULL) 			free(appCF);
    if(PiDi != NULL) 			free(PiDi);
	if(rmsd_matrix != NULL) 	free(rmsd_matrix);
	if(density_matrix != NULL) 	free(density_matrix);
	if(nearest_center != NULL)	free(nearest_center);
	if(distance_matrix != NULL) free(distance_matrix);
    if(cartesian_coord != NULL) free(cartesian_coord);
	if(assigned_cluster != NULL)free(assigned_cluster);
}


void assign_cluster_from_density_neighborhood(int i, int* nUnclustered, int* nOutliers, int* density_matrix, float* distance_matrix, int* nearest_center, int* assigned_cluster)
{
    // rewrite this function to assign cluster recursively
}


void QuickSort_PiDi(chromosome* chrom, double* appCF, float* distance_matrix, int* density_matrix, int* assigned_cluster, int* nearest_center, float* PiDi, float* coord, float* rmsd, int num_chrom, int beg, int end)
{
	int l,r,p;
	float pivot;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		pivot = PiDi[p];

		while(1)
		{
			while( (l<=r) && QS_DSC(PiDi[l],pivot) <= 0 ) ++l;
			while( (l<=r) && QS_DSC(PiDi[r],pivot)  > 0 ) --r;
			
			if (l > r) break;

			swap_elements(	&chrom[l],&appCF[l],&distance_matrix[l],&density_matrix[l], &assigned_cluster[l], &nearest_center[l], &PiDi[l], &coord[l*3*MAX_ATM_HET], &rmsd[num_chrom*l],
							&chrom[r],&appCF[r],&distance_matrix[r],&density_matrix[r], &assigned_cluster[r], &nearest_center[r], &PiDi[r], &coord[r*3*MAX_ATM_HET], &rmsd[num_chrom*r], num_chrom);
			if (p == r) p=l;
			++l;--r;
		}
		swap_elements(	&chrom[p],&appCF[p],&distance_matrix[p],&density_matrix[p], &assigned_cluster[p], &nearest_center[p], &PiDi[p], &coord[p*3*MAX_ATM_HET], &rmsd[p*num_chrom],
						&chrom[r],&appCF[r],&distance_matrix[r],&density_matrix[r], &assigned_cluster[r], &nearest_center[r], &PiDi[r], &coord[r*3*MAX_ATM_HET], &rmsd[r*num_chrom], num_chrom);
		
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_PiDi(chrom, appCF, distance_matrix, density_matrix, assigned_cluster, nearest_center, PiDi, coord, rmsd, num_chrom, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_PiDi(chrom, appCF, distance_matrix, density_matrix, assigned_cluster, nearest_center, PiDi, coord, rmsd, num_chrom, l, end);
			end = r;
		}
	}
}

// this function is used to reorder cluster by decreasing density
void QuickSort_Density(chromosome* chrom, double* appCF, float* distance_matrix, int* density_matrix, int* assigned_cluster, int* nearest_center, float* PiDi, float* coord, float* rmsd, int num_chrom, int beg, int end)
{
	int l,r,p;
	int pivot;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		pivot = density_matrix[p];

		while(1)
		{
			while( (l<=r) && QS_DSC(density_matrix[l],pivot) <= 0 ) ++l;
			while( (l<=r) && QS_DSC(density_matrix[r],pivot)  > 0 ) --r;
			
			if (l > r) break;

			swap_elements(	&chrom[l],&appCF[l],&distance_matrix[l],&density_matrix[l], &assigned_cluster[l], &nearest_center[l], &PiDi[l], &coord[l*3*MAX_ATM_HET], &rmsd[num_chrom*l],
							&chrom[r],&appCF[r],&distance_matrix[r],&density_matrix[r], &assigned_cluster[r], &nearest_center[r], &PiDi[r], &coord[r*3*MAX_ATM_HET], &rmsd[num_chrom*r], num_chrom);
			if (p == r) p=l;
			++l;--r;
		}
		swap_elements(	&chrom[p],&appCF[p],&distance_matrix[p],&density_matrix[p], &assigned_cluster[p], &nearest_center[p], &PiDi[p], &coord[p*3*MAX_ATM_HET], &rmsd[p*num_chrom],
						&chrom[r],&appCF[r],&distance_matrix[r],&density_matrix[r], &assigned_cluster[r], &nearest_center[r], &PiDi[r], &coord[r*3*MAX_ATM_HET], &rmsd[r*num_chrom], num_chrom);
		--r;

		if( (r-beg) < (end-l) )
		{
			QuickSort_Density(chrom, appCF, distance_matrix, density_matrix, assigned_cluster, nearest_center, PiDi, coord, rmsd, num_chrom, beg, r);
			beg = l;
		}
		else
		{
			QuickSort_Density(chrom, appCF, distance_matrix, density_matrix, assigned_cluster, nearest_center, PiDi, coord, rmsd, num_chrom, l, end);
			end = r;
		}
	}
}


void swap_elements(chromosome* CHROMx, double* appCFx, float* DISTx, int* DENSITYx, int* CLUSx, int* CENx, float* PiDix, float* COORx, float* RMSDx, chromosome* CHROMy, double* appCFy, float* DISTy, int* DENSITYy, int* CLUSy, int* CENy, float* PiDiy, float* COORy, float* RMSDy, int num_chrom)
{
	// temporary variables 
	int i,j;
	double appCFt = 0.0;
	chromosome CHROMt;
	float DISTt = 0.0f, PiDit = 0.0f;
	int DENSITYt = 0, CLUSt = 0, CENt = -1;
	float* COORt = (float*) malloc(MAX_ATM_HET*3*sizeof(float));
	float* RMSDt  = (float*) malloc(num_chrom*sizeof(float));
	//mem check
	if(!COORt || !RMSDt) {fprintf(stderr,"ERROR: The Partition Function is NULL in the clustering step.\n"); Terminate(2);}
	//mem init
	memset(COORt, 0.0f, MAX_ATM_HET*3);
	memset(RMSDt, 0.0f, num_chrom);
	//mem swap
	appCFt = *appCFx; *appCFx = *appCFy; *appCFy = appCFt;
	CHROMt = *CHROMx; *CHROMx = *CHROMy; *CHROMy = CHROMt;
	DISTt = *DISTx; *DISTx = *DISTy; *DISTy = DISTt;
	DENSITYt = *DENSITYx; *DENSITYx = *DENSITYy; *DENSITYy = DENSITYt;
	CLUSt = *CLUSx; *CLUSx = *CLUSy; *CLUSy = CLUSt;
	CENt = *CENx; *CENx = *CENy; *CENy = CENt;
	PiDit = *PiDix; *PiDix = *PiDiy; *PiDiy = PiDit;
	for(i = 0; i < num_chrom; ++i)		{RMSDt[i] = RMSDx[i]; RMSDx[i] = RMSDy[i]; RMSDy[i] = RMSDt[i];}
	for(i = 0; i < 3*MAX_ATM_HET; ++i)	{COORt[i] = COORx[i]; COORx[i] = COORy[i]; COORy[i] = COORt[i];}
    if(COORt) free(COORt);
    if(RMSDt) free(RMSDt);
}

float calculate_stddev(float* PiDi, int num_chrom)
{
	int i;
	float mean = calculate_mean(PiDi, num_chrom);
	float sqrtot = 0.0f;
    
    // printf("mean PiDi = %g\n", mean);
	for(i = 0; i < num_chrom; ++i) sqrtot += pow( (mean - PiDi[i]),2.0 );
	return sqrtf(sqrtot/(float)num_chrom);
}

float calculate_mean(float* PiDi, int num_chrom)
{
	int i;
	float mean = 0.0f;
	// printf("float max = %g\n", FLT_MAX);
    for(i = 0; i < num_chrom; ++i) mean += PiDi[i];
	return mean /= (float)num_chrom;
}
