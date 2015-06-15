#ifndef FOPTICS_H
#define FOPTICS_H

#include "gaboom.h"
#include "boinc.h"

#define UNDEFINED_DIST -0.1f // Defined in FOPTICS as > than +INF
#define minPts 5
/*****************************************\
			RandomProjections
\*****************************************/
class RandomProjectedNeighborsAndDensities
{
	public:
		RandomProjectedNeighborsAndDensities(); // Constructor (publicly called from FlexAID *_cluster.cxx)
		void computeSetBounds();
	
	private:
		int N;
		int nDimensions;
		static const int RandomSeed;
		static const int logOProjectionConstant = 20;
		static const double sizeTolerance = (double) 2.0/3.0;
		std::vector< chromosome* > chroms;
		int minSplitSize;

	protected:
		int nProject1D;
		std::vector< std::vector< int > > splitsets;
		std::vector< std::vector< double > > projectedPoints;
};

/*****************************************\
				FastOPTICS
\*****************************************/
class FastOPTICS
{
	friend class RandomProjectedNeighborsAndDensities;
	
	public:
		FastOPTICS(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom); // Constructor (publicly called from FlexAID's *_cluster.cxx)
	
	private:
		// FlexAID specific attributes
		int N;	// N : number of chromosomes to cluster
		FA_Global* FA;			// pointer to FA_Global struct
		GB_Gblobal* GB;			// pointer to GB_Global struct
		VC_Global* VC;			// pointer to VC_Global struct
		chromosome* chrom;		// pointer to chromosomes' array
		genlim* gene_lim;		// pointer to gene_lim genlim array (useful for bondaries defined for each gene)
		atom* atoms;			// pointer to atoms' array
		resid* residue;			// pointer to residues' array
		gridpoint* cleftgrid;	// pointer to gridpoints' array (defining the total search space of the simulation)
		
		// FOPTICS algorithm attributes
		static int iOrder;
		std::vector< int > order;
		std::vector< double > reachDist;
		std::vector< bool > processed;
		std::vector< double > inverseDensities;
		std::vector< chromosome* > points;
		std::vector< std::vector< chromosome* > > neighs;
		
		// private methods
		void Initialize(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom); // Initialize FastOPTICS private attributes from FlexAID structs
};

#endif