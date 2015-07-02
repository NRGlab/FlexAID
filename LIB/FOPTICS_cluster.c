#include "FOPTICS.h"

void FOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom)
{
	// create Cluster or BindingMode vetor to send to FastOPTICS constructor
	BindingPopulation::BindingPopulation Population(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
	// calling FastOPTICS constructor
	// 	1. Partition Sets using Random Vectorial Projections
	// 	2. Calculate Neighborhood
	// 	3. Calculate reachability distance
	// 	4. Compute the Ordering of Points To Identify Cluster Structure (OPTICS)
	// 	5. Populate BindingPopulation::Population after analyzing OPTICS
	FastOPTICS::FastOPTICS(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom, Population&);

}
