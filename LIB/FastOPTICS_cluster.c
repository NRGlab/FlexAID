#include "FOPTICS.h"

void FastOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
   // int minPoints = static_cast<int>( floor(nChrom * 0.0005) );
    int minPoints = 8;
    // (minPoints < 3*FA->num_het_atm) ? minPoints = minPoints : minPoints = 3*FA->num_het_atm;
	
    // BindingPopulation() : BindingPopulation constructor *non-overridable*
    BindingPopulation::BindingPopulation Population1(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    BindingPopulation::BindingPopulation Population2(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    BindingPopulation::BindingPopulation Population3(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    BindingPopulation::BindingPopulation Population4(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
	BindingPopulation::BindingPopulation Population5(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    
    // FastOPTICS() : calling FastOPTICS constructors
    FastOPTICS::FastOPTICS Algo1(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population1, minPoints);
    minPoints = std::floor(minPoints * 1.5);
    FastOPTICS::FastOPTICS Algo2(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population2, minPoints);
    minPoints = std::floor(minPoints * 1.5);
    FastOPTICS::FastOPTICS Algo3(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population3, minPoints);
    minPoints = std::floor(minPoints * 1.5);
    FastOPTICS::FastOPTICS Algo4(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population4, minPoints);
    minPoints = std::floor(minPoints * 1.5);
    FastOPTICS::FastOPTICS Algo5(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population5, minPoints);
    
    // 	1. Partition Sets using Random Vectorial Projections
    // 	2. Calculate Neighborhood
    // 	3. Calculate reachability distance
    // 	4. Compute the Ordering of Points To Identify Cluster Structure (OPTICS)
    // 	5. Populate BindingPopulation::Population after analyzing OPTICS
    Algo1.Execute_FastOPTICS(end_strfile, tmp_end_strfile);
    Algo2.Execute_FastOPTICS(end_strfile, tmp_end_strfile);
    Algo3.Execute_FastOPTICS(end_strfile, tmp_end_strfile);
    Algo4.Execute_FastOPTICS(end_strfile, tmp_end_strfile);
    Algo5.Execute_FastOPTICS(end_strfile, tmp_end_strfile);

    Algo1.output_OPTICS(end_strfile, tmp_end_strfile);
    Algo2.output_OPTICS(end_strfile, tmp_end_strfile);
    Algo3.output_OPTICS(end_strfile, tmp_end_strfile);
    Algo4.output_OPTICS(end_strfile, tmp_end_strfile);
    Algo5.output_OPTICS(end_strfile, tmp_end_strfile);
    
    std::cout << "Size of Population 1 is " << Population1.get_Population_size() << " Binding Modes." << std::endl;
    std::cout << "Size of Population 2 is " << Population2.get_Population_size() << " Binding Modes." << std::endl;
    std::cout << "Size of Population 3 is " << Population3.get_Population_size() << " Binding Modes." << std::endl;
    std::cout << "Size of Population 4 is " << Population4.get_Population_size() << " Binding Modes." << std::endl;
    std::cout << "Size of Population 5 is " << Population5.get_Population_size() << " Binding Modes." << std::endl;
    
    // output FA->max_result BindingModes
    Population1.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo1.get_minPoints());
    Population2.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo2.get_minPoints());
    Population3.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo3.get_minPoints());
    Population4.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo4.get_minPoints());
    Population5.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo5.get_minPoints());
    printf("-- end of FastOPTICS_cluster --\n");
}