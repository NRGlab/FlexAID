#ifndef FOPTICS_H
#define FOPTICS_H

#include "gaboom.h"
#include "boinc.h"

#define UNDEFINED_DIST -0.1f // Defined in FOPTICS as > than +INF

/*****************************************\
			RandomProjections
\*****************************************/
class RandomProjectedNeighborsAndDensities
{
	public:
		RandomProjectedNeighborsAndDensities(); // Constructor (publicly called from FlexAID *_cluster.cxx)

	private:
		int N;
		int nDimensions;
		static const int RandomSeed;
		static const int logOProjectionConstant = 20;
		const static  double sizeTolerance = (double) 2.0/3.0;
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
		FastOPTICS(); // Constructor (publicly called from FlexAID's *_cluster.cxx)
		Initialize(); // Initialize FastOPTICS private attributes from FlexAID structs
	private:
		// FlexAID specific attributes
		int N;					// N : number of chromosomes to cluster
		FA_Global* FA;			// pointer to FA_Global struct
		GB_Gblobal* GB;			// pointer to GB_Global struct
		VC_Global* VC;			// pointer to VC_Global struct
		chromosome* chrom;		// pointer to chromosomes' array
		genlim* gene_lim;		// pointer to gene_lim genlim array (useful for bondaries defined for each gene)
		atom* atoms;			// pointer to atoms' array
		resid* residue;			// pointer to residues' array
		gridpoint* cleftgrid;	// pointer to gridpoints' array (defining the total search space of the simulation)
		
		// FOPTICS algorithm attributes
		int iOrder;
		std::vector< int > order;
		std::vector< double > reachDist;
		std::vector< bool > processed;
		std::vector< double > inverseDensities;
		std::vector< chromosome* > chroms;
		std::vector< std::vector< chromosome* > > neighs;
};

#endif