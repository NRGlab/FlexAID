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
    
    std::cout << "Size of Population is " << Population.get_Population_size() << " Binding Modes." << std::endl;
    
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

void entropy_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp)
{
    // minPoints will be used for the call to ColonyEnergy
    // "       " will also be used for output_Population as a suffix (TO UNIFORMIZE)
    int minPoints = 15;
    int nClusters = 0;
    std::clock_t start, check1, check2, check3, check4;
    double duration = 0.0;
    
    start = std::clock();

    BindingPopulation Population = BindingPopulation(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom);

    check1 = std::clock();
    duration = ( check1 - start ) / (double) CLOCKS_PER_SEC;
    std::cout<<"BindingPopulation() : "<< duration <<std::endl;
    
    // call to ColonyEnergy class will serve to get neighbors for each chrom
    // this will be useful later to build BindingModes in BindingPopulation
    // ColonyEnergy Algo = ColonyEnergy(FA, GB, VC, chrom, gene_lim, atoms, residue, cleftgrid, nChrom, Population, minPoints);
    // Algo.Execute_ColonyEnergy(end_strfile, tmp_end_strfile);
    
    check2 = std::clock();
    duration = ( check2 - check1 ) / (double) CLOCKS_PER_SEC;
    std::cout<<"ColonyEnergy : "<< duration <<std::endl;
    
    // iterates over all Poses to run the CF clustering algo
    for(std::vector<Pose>::iterator iPose = Population.Poses.begin(); iPose != Population.Poses.end(); ++iPose)
    {
        if(!iPose->processed) // pose is still unclustered
        {
            BindingMode mode(&Population);
            iPose->order = nClusters;
            iPose->processed = true;
            mode.add_Pose(*iPose);

            for(std::vector<Pose>::iterator jPose = iPose+1; jPose != Population.Poses.end(); ++jPose)
            {   
                if(!jPose->processed) // pose is still unclustered
                {
                    float rmsd = Population.compute_distance(*iPose, *jPose);
                    // float rmsd2 = calc_rmsd_chrom(FA,GB,chrom,gene_lim,atoms,residue,cleftgrid,GB->num_genes,iPose->chrom_index,jPose->chrom_index, NULL, NULL, true);
                    
                    if( rmsd <= FA->cluster_rmsd )
                    {
                        jPose->order = nClusters;
                        jPose->reachDist = rmsd;
                        jPose->processed = true;
                        mode.add_Pose(*jPose);
                    }
                }
            }
            
            check3 = std::clock();
            duration = ( check3 - check2 ) / (double) CLOCKS_PER_SEC;
            std::cout<<"CF cluster : "<< duration <<std::endl;
            
            // neighbors-adding section
            // std::vector<int> neighs = Algo.get_neighbors_for_chrom(iPose->chrom_index);
            
            // iterating over neighbors to add them to BindingMode
            // for(std::vector<int>::iterator it = neighs.begin(); it != neighs.end(); ++it)
            // {
            //     // continue to next iter if the chrom_index of the Pose is already in BindingMode
            //     if(mode.isPoseInBindingMode(*it)) continue;
            //     // need to find the appropriate Pose by index (maybe all chroms were not built as Poses in 1.)
            //     for(std::vector<Pose>::iterator itPose = Population.Poses.begin(); itPose != Population.Poses.end(); itPose++)
            //     {
            //         if( *it == itPose->chrom_index )
            //         {
            //             itPose->order = nClusters;
            //             itPose->processed = true;
            //             // itPose->reachDist = Population.compute_distance(*iPose, *itPose);
            //             /*if( mode.isPoseAggregable(*itPose) )*/ mode.add_Pose(*itPose);
            //             break;
            //         }
            //     }
            // }
            
            check4 = std::clock();
            duration = ( check4 - check3 ) / (double) CLOCKS_PER_SEC;
            std::cout<<"Neighbors : "<< duration <<std::endl;
            
            Population.add_BindingMode(mode);
            nClusters++;
        }

        if(nClusters == FA->max_results) { break; } // stop clustering when enough results are generated
    }
    // Population.Classify_BindingModes();
    Population.output_Population(FA->max_results, end_strfile, tmp_end_strfile, dockinp, gainp, minPoints); // minPoints = 0 to call the same function without overloading/modifying output_Population 
}
