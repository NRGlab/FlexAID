#include "FOPTICS.h"

void FastOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
	// BindingPopulation() : BindingPopulation constructor *non-overridable*
	BindingPopulation::BindingPopulation Population(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
	
    // FastOPTICS() : calling FastOPTICS constructor
    FastOPTICS::FastOPTICS Algo(FA, GB, VC, chrom, gene_lim, cleftgrid, nChrom, Population);
    
    // 	1. Partition Sets using Random Vectorial Projections
    // 	2. Calculate Neighborhood
    // 	3. Calculate reachability distance
    // 	4. Compute the Ordering of Points To Identify Cluster Structure (OPTICS)
    // 	5. Populate BindingPopulation::Population after analyzing OPTICS
	Algo.Execute_FastOPTICS();
    Algo.output_OPTICS(end_strfile, tmp_end_strfile);
    std::cout << "Size of Population is " << Population.get_Population_size() << " Binding Modes." << std::endl;
    
    // output FA->max_result BindingModes
    Population.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp);
//    Population.output_Population(25, end_strfile, tmp_end_strfile, dockinp, gainp);
}