#include "gaboom.h"
#include "boinc.h"

// the following macro is used to retrieve the appropriate value from a matrix 
#define I(i, j, dim) ( (i<=j) ? (i*dim+j) : (j*dim+i) )
#define V(mat, i, j, dim)  ( if(mat) { (i <= j) ? mat[i*dim+j] : mat[j*dim+i] } )

void density_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int memchrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
	// The following steps will be executed prior to assign clusters
	/*	(1) Build RSMD matrix between each chrom [dij]
		(2) Build Local Density matrix of each chrom (Eq. 1) [pij] : pi = ∑jCHI(dij - d_cutoff)
		(3) Build Distance from Density Centers matrix (Eq. 2) [∂i] = min(dij) for pj > pi
		(4) Order chrom in decreasing ∂i : ∂i is much larger for chrom which are local/global maxima in density
		(5) Cluster centers are recognized as chrom with anomalously large ∂i values
	*/

	// variables initialization
	cfstr cf;										/* Complementarity Function value */
	int i = 0, j = 0;								/* indexing iterators */
	int dim = ( memchrom * (memchrom + 1) ) / 2;	/* 1D array dimension needed to emulate storage of a symmetric 2D matrix */
	float* rmsd_matrix = NULL;						/*  (1) */
	int* density_matrix = NULL;						/*  (2) */ 
	float* distance_matrix = NULL;					/*  (3) */ 
	bool Hungarian = false; 						/* Bool indicating whether to use symmetry corrected RMSD or not */

	// memory allocation
	rsmd_matrix = (float*)malloc(dim*sizeof(float));
	density_matrix = (int*)malloc(dim*sizeof(int));
	distance_matrix = (float*)malloc(dim*sizeof(float));

	// memory check
	if(!rmsd_matrix || !density_matrix || !distance_matrix)
	{
		fprintf(stderr,"ERROR: memory allocation error for clusters\n");
		Terminate(2);
	}

	// matrices initialization
	for(i = 0; i <= dim; ++i)
	{
		density_matrix[i] = 0;
		rmsd_matrix[i] = 0.0f;
		distance_matrix[i] = 0.0f;
	}

	// (1) build RMSD matrix between each chrom
	for(i = 0; i <= memchrom; ++i)
	{
		for()
	}
}