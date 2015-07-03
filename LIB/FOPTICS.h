#ifndef FOPTICS_H
#define FOPTICS_H

#include "BindingMode.h"
#include <utility>
#include <queue>

#define UNDEFINED_DIST -0.1f // Defined in FOPTICS as > than +INF
#define isUndefinedDist(a) ((a - UNDEFINED_DIST) <= FLT_EPSILON)

struct ClusterOrdering
{
	int objectID;
	int predecessorID;
	float reachability;

	ClusterOrdering(int, int, float);

	inline bool const operator==(const ClusterOrdering& rhs);
	inline bool const operator< (const ClusterOrdering& rhs);

};

struct ClusterOrderingComparator
{
	inline bool operator() ( const ClusterOrdering& pose1, const ClusterOrdering& pose2 )
	{
		if(pose1.reachability > pose2.reachability)
			return 1;
		else if(pose1.reachability < pose2.reachability)
			return 1;
		if(pose1.objectID > pose2.objectID)
			return -1;
		else if(pose1.objectID < pose2.objectID)
			return 1;
		// if nothing else is true, return 0
		return 0;
	}
};


/*****************************************\
				FastOPTICS
\*****************************************/
class FastOPTICS
{
	friend class RandomProjectedNeighborsAndDensities;
	
	public:
		FastOPTICS(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, BindingPopulation&); // Constructor (publicly called from FlexAID's *_cluster.cxx)
	
	private:
		// FlexAID specific attributes
		int N;			// N : number of chromosomes to cluster
		int minPoints;	// minPts : minimal number of neighbors (only parameter in FOPTICS)
		int nDimensions;
		
		// FOPTICS algorithm attributes
		static int iOrder;
		std::vector< int > order;
		std::vector< float > reachDist;
		std::vector< bool > processed;
		std::vector< float > inverseDensities;
		std::vector< std::pair<chromosome*,std::vector<float> > > points;
		std::vector< std::vector< int > > neighbors;

		// BindingPopulation is used for clustering purposed
		BindingPopulation* Population;
		
		// private methods
		void Initialize(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom); // Initialize FastOPTICS private attributes from FlexAID structs
		void ExpandClusterOrder(int);
		// void Clusterize();
	

	protected:
		// protected methods to be used by RandomProjectedNeighborsAndDensities::methods()
		static FA_Global* FA;			// pointer to FA_Global struct
		static GB_Global* GB;			// pointer to GB_Global struct
		static VC_Global* VC;			// pointer to VC_Global struct
		static chromosome* chroms;		// pointer to chromosomes' array
		static genlim* gene_lim;		// pointer to gene_lim genlim array (useful for bondaries defined for each gene)
		static atom* atoms;				// pointer to atoms' array
		static resid* residue;			// pointer to residues' array
		static gridpoint* cleftgrid;	// pointer to gridpoints' array (defining the total search space of the simulation)
		std::vector<float> 				Vectorized_Chromosome(chromosome* chrom);
		float							compute_distance(std::pair< chromosome*,std::vector<float> > &, std::pair< chromosome*,std::vector<float> > &);
};

/*****************************************\
			RandomProjections
\*****************************************/
class RandomProjectedNeighborsAndDensities
{
	friend class FastOPTICS;
	
	public:
		RandomProjectedNeighborsAndDensities(std::vector< std::pair< chromosome*,std::vector<float> > >&, int, FastOPTICS*); // Constructor (publicly called from FlexAID *_cluster.cxx)
	
	private:
		// provate attributes
		FastOPTICS* top;
		int N;
		int nDimensions;
		int minSplitSize;
		static const int logOProjectionConstant = 20;
		static float sizeTolerance;
		std::vector< std::pair<chromosome*,std::vector<float> > > points;

		// private methods
		void 								SplitUpNoSort(std::vector<int>&, int);
		void 								quicksort_concurrent_Vectors(std::vector<float>& data, std::vector<int>& index, int begin, int end);
		void 								swap_element_in_vectors(std::vector<float>::iterator, std::vector<float>::iterator, std::vector<int>::iterator , std::vector<int>::iterator);
		
	protected:
		// protected attributes (accessible via FastOPTICS class)
		int nProject1D;				// CONSTANT c0
		int nPointsSetSplits;		// CONSTANT c1
		std::vector< std::vector< int > > 	splitsets;
		std::vector< std::vector< float > > projectedPoints;
		// protected methods (accessible via FastOPTICS class)
		void 								computeSetBounds(std::vector< int >&);
		std::vector<float> 					getInverseDensities();
		std::vector< std::vector< int > > 	getNeighbors();
		std::vector<float> 					Randomized_Normalized_Vector();
		int 								Dice();
};
#endif