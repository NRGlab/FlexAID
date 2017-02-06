#include "FOPTICS.h"
#include "ColonyEnergy.h"

/*
    The following functions will be called from top.c 
 */
void FastOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
    int minPoints = 15;
    // (minPoints < 3*FA->num_het_atm) ? minPoints = minPoints : minPoints = 3*FA->num_het_atm;
	
    // BindingPopulation() : BindingPopulation constructor *non-overridable*
    // BindingPopulation::BindingPopulation Population(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    BindingPopulation Population = BindingPopulation(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);
    
    // FastOPTICS() : calling FastOPTICS constructors
    FastOPTICS Algo = FastOPTICS(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population, minPoints);
    
    // 	1. Partition Sets using Random Vectorial Projections
    // 	2. Calculate Neighborhood
    // 	3. Calculate reachability distance
    // 	4. Compute the Ordering of Points To Identify Cluster Structure (OPTICS)
    // 	5. Populate BindingPopulation::Population after analyzing OPTICS
    Algo.Execute_FastOPTICS(end_strfile, tmp_end_strfile);

    // Algo.output_OPTICS(end_strfile, tmp_end_strfile);

    // output the 3D poses ordered with Fast OPTICS (done only once for the purpose as the order should not change)
    // Algo.output_3d_OPTICS_ordering(end_strfile, tmp_end_strfile);
    
    std::cout << "Size of Population 1 is " << Population.get_Population_size() << " Binding Modes." << std::endl;
    
    // output FA->max_result BindingModes
    Population.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo.get_minPoints());
    // printf("-- end of FastOPTICS_cluster --\n");
}


void ColonyEnergy_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
    int minPoints = 15;

    BindingPopulation Population = BindingPopulation(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);

    ColonyEnergy Algo = ColonyEnergy(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population, minPoints);

    Algo.Execute_ColonyEnergy(end_strfile, tmp_end_strfile);

    Population.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, Algo.get_minPoints());
}

void Entropy_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
    // minPoints will be used for the call to ColonyEnergy
    // "       " will also be used for output_Population as a suffix (TO UNIFORMIZE)
    int minPoints = 15;

    BindingPopulation Population = BindingPopulation(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftfrid, nChrom);

    // call to ColonyEnergy class is still incertain for the current implementation
    ColonyEnergy Algo = ColonyEnergy(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population, minPoints);
    
    // multiples options to consider :
    // 1. Aggregate by CF then call ColonyEnergy (need to modify ColonyEnergy)
    // 2. Process ColonyEnergy
    std::vector<Pose> poses;
    for(int i = 0; i < nChrom; ++i)
    {   
        BindingMode mode(&Population);

        if(chrom[i].app_evalue < 0)
        {
            std::vector<float> vPose(Algo.Vectorized_Cartesian_Coordinates(i));
            Pose pose(&chrom[i], i, -1, 0.0f, FA->temperature, vPose);
            poses.push_back(pose);
        }
    } 

    Population.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, minPoints); // minPoints = 0 to call the same function without overloading/modifying output_Population 
}