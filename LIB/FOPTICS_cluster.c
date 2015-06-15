#include "FOPTICS.h"

void FOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom)
{
	// create Cluster or BindingMode vetor to send to FastOPTICS constructor
	
	// calling FastOPTICS constructor
	FastOPTICS::FastOPTICS(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
}