#ifndef COLONYENERGY_H
#define COLONYENERGY_H

#include "BindingMode.h"
#include <utility>
#include <iostream>
#include <sstream> // for ostringstream
#include <string>
#include <fstream>
#include <queue>
#include <cmath>


/*****************************************\
				FastOPTICS
\*****************************************/
class ColonyEnergy
{
	// multiple random projections partitioning algorithm-related class
	friend class RandomProjectedNeighborsColonyEnergy;
	// entropy-related classes to integrate with FlexAID
	friend class BindingMode;
    friend class BindingPopulation;
	
	public:
		// Constructor (publicly called from FlexAID's *_cluster.cxx)
		explicit 	ColonyEnergy(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, BindingPopulation&, int nPoints);
		// Main FastOPTICS function to execute the algorithm
		void 		Execute_ColonyEnergy(char* end_strfile, char* tmp_end_strfile);
		
		// population classification methods (linked to BindingPopulation* Population)


        // distance computing methods
		float 				compute_distance(std::pair< chromosome*,std::vector<float> > &, std::pair< chromosome*,std::vector<float> > &);
		float 				compute_vect_distance(std::vector<float> a, std::vector<float> b);
		// getters methods
		int 				get_minPoints();
		std::vector<int> 	get_neighbors_for_chrom(int);
		std::vector<float> 	Vectorized_Chromosome(chromosome* chrom);
		std::vector<float>	Vectorized_Cartesian_Coordinates(int chrom_index);

	private:
		// FlexAID specific attributes
		int 	N;				// N : number of chromosomes to cluster
		int 	minPoints;		// minPts : minimal number of neighbors (only parameter in FOPTICS)
		int 	nDimensions;	// number of dimensions of the vectorized Pose
		float 	dist_threshold; // contains the value of this->FA->cluster_rmsd*(2 - RandomProjectedNeighborsColonyEnergy::sizeTolerance)
		
		// FOPTICS algorithm attributes
		// std::vector< int > order;
		// std::vector< float > reachDist;
		std::vector< bool > processed;
		std::vector< float > inverseDensities;
		std::vector< std::pair<chromosome*,std::vector<float> > > points;
		std::vector< std::vector< int > > neighbors;
        
        std::vector< Pose > OPTICS;
        
		// BindingPopulation is used for clustering purposed
		BindingPopulation* Population;
		
		// private methods
		void 	ExpandClusterOrder(int);
		void 	normalizeDistances();
	
	
	protected:
		// protected methods to be used by RandomProjectedNeighborsColonyEnergy::methods()
			  	 FA_Global* FA;		// pointer to FA_Global struct
		/*const*/GB_Global* GB;		// pointer to GB_Global struct
		/*const*/VC_Global* VC;		// pointer to VC_Global struct
		/*const*/chromosome* chroms;	// pointer to chromosomes' array
		/*const*/genlim* gene_lim;		// pointer to gene_lim genlim array (useful for bondaries defined for each gene)
		atom* atoms;				// pointer to atoms' array
		resid* residue;				// pointer to residues' array
		/*const*/gridpoint* cleftgrid;	// pointer to gridpoints' array (defining the total search space of the simulation)
};

/*****************************************\
			RandomProjections
\*****************************************/
class RandomProjectedNeighborsColonyEnergy
{
	friend class ColonyEnergy;
	
	public:
		explicit RandomProjectedNeighborsColonyEnergy(std::vector< std::pair< chromosome*,std::vector<float> > >&, int, ColonyEnergy*); // Constructor (publicly called from FlexAID *_cluster.cxx)
	
	private:
		// private attributes
		ColonyEnergy* top;
		int N;
		int nDimensions;
		int minSplitSize;
		static const int logOProjectionConstant;
		static float sizeTolerance;
		std::vector< std::pair<chromosome*,std::vector<float> > > points;

		// private methods
		void 								SplitUpNoSort(std::vector<int>&, int);
		void 								quicksort_concurrent_Vectors(std::vector<float>& data, std::vector<int>& index, int begin, int end);
		void 								swap_element_in_vectors(std::vector<float>::iterator, std::vector<float>::iterator, std::vector<int>::iterator , std::vector<int>::iterator);
		
	protected:
		// protected attributes (accessible via ColonyEnergy class)
		int nProject1D;				// CONSTANT c0
		int nPointsSetSplits;		// CONSTANT c1
		std::vector< std::vector< int > > 	splitsets;
		std::vector< std::vector< float > > projectedPoints;
		// protected methods (accessible via ColonyEnergy class)
		void 								computeSetBounds(std::vector< int >&);
		void								getInverseDensities(std::vector<float> &);
		void								getNeighbors(std::vector< std::vector< int > > &);
		std::vector<float> 					Randomized_InternalCoord_Vector();
        std::vector<float>                  Randomized_CartesianCoord_Vector();
        std::vector<float>                  Randomly_Selected_Chromosome();
};
#endif