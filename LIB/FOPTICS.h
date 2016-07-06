#ifndef FOPTICS_H
#define FOPTICS_H

#include "BindingMode.h"
#include <utility>
#include <iostream>
#include <sstream> // for ostringstream
#include <string>
#include <fstream>
#include <queue>
#include <cmath>



// random die &&  Float comparators
int roll_die();
bool definitelyGreaterThan(float a, float b, float epsilon);
bool definitelyLessThan(float a, float b, float epsilon);

float normalize_IC_interval(const genlim* gene_lim, float dist);

struct ClusterOrdering
{
	int objectID;
	int predecessorID;
	float reachability;

	ClusterOrdering(int, int, float);

	inline bool const operator==(const ClusterOrdering& rhs);
	inline bool const operator< (const ClusterOrdering& rhs);
	inline bool const operator> (const ClusterOrdering& rhs);

};

struct ClusterOrderingComparator
{
	// inline int operator() ( const ClusterOrdering& pose1, const ClusterOrdering& pose2 )
	inline bool operator() ( const ClusterOrdering& pose1, const ClusterOrdering& pose2 )
	{
		if(pose1.reachability > pose2.reachability || isUndefinedDist(pose1.reachability))
			return true;
		else if(pose1.reachability < pose2.reachability)
			return false;

		// if(pose1.predecessorID > pose2.predecessorID)
		// 	return true;
		// else if(pose1.predecessorID < pose2.predecessorID)
		// 	return false;

		if(pose1.objectID > pose2.objectID)
			return false;
		else if(pose1.objectID < pose2.objectID)
			return true;

		// if nothing else is true, return 0
		return false;
	}
};

// priority_queue container alllowing iterations throught the priority_queue
template <class T, class S, class C>S& Container(priority_queue<T, S, C>& q)
{
    struct IterableQueue : private priority_queue<T, S, C>
    {
        static S& Container(priority_queue<T, S, C>& q)
        {
            return q.*&IterableQueue::c;
        }
    };
    return IterableQueue::Container(q);
}

/*****************************************\
				FastOPTICS
\*****************************************/
class FastOPTICS
{
	// multiple random projections partitioning algorithm
	friend class RandomProjectedNeighborsAndDensities;
	// entropy-related classes to integrate with FlexAID
	friend class BindingMode;
    friend class BindingPopulation;
	
	public:
		// Constructor (publicly called from FlexAID's *_cluster.cxx)
		explicit 	FastOPTICS(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, BindingPopulation&, int nPoints);
		// Main FastOPTICS function to execute the algorithm
		void 		Execute_FastOPTICS(char* end_strfile, char* tmp_end_strfile);
		// void		update_ClusterOrdering_PriorityQueue_elements(int, std::priority_queue< ClusterOrdering, std::vector<ClusterOrdering>, ClusterOrderingComparator::ClusterOrderingComparator > &);
		// population classification methods (linked to BindingPopulation* Population)
				// 	classify the Pose into the appropriate BindingMode in the BindingPopulation
		bool 		Classify_Pose(Pose& pose);
				// 	classify the whole BindingPopulation into BindingModes (this method might call Classify_Pose())
		void 		Classify_Population();
		// output-related methods
        void 		output_OPTICS(char* end_strfile, char* tmp_end_strfile);
        void 		output_3d_OPTICS_ordering(char* end_strfile, char* tmp_end_strfile);
        // distance computing methods
		float 		compute_distance(std::pair< chromosome*,std::vector<float> > &, std::pair< chromosome*,std::vector<float> > &);
		float 		compute_vect_distance(std::vector<float> a, std::vector<float> b);
		// getters methods
		int 		get_minPoints();

	private:
		// FlexAID specific attributes
		int 	N;			// N : number of chromosomes to cluster
		int 	minPoints;	// minPts : minimal number of neighbors (only parameter in FOPTICS)
		int 	nDimensions;
		float 	dist_threshold;
		// FOPTICS algorithm attributes
		int iOrder;
		std::vector< int > order;
		std::vector< float > reachDist;
		std::vector< bool > processed;
		std::vector< float > inverseDensities;
		std::vector< std::pair<chromosome*,std::vector<float> > > points;
		std::vector< std::vector< int > > neighbors;
        
        std::vector< Pose > OPTICS;
        // std::priority_queue< Pose, std::vector<Pose>, PoseClassifier::PoseClassifier > OPTICS;
		// BindingPopulation is used for clustering purposed
		BindingPopulation* Population;
		
		// private methods
		void 	ExpandClusterOrder(int);
		void 	normalizeDistances();
	
	
	protected:
		// protected methods to be used by RandomProjectedNeighborsAndDensities::methods()
			  	 FA_Global* FA;		// pointer to FA_Global struct
		/*const*/GB_Global* GB;		// pointer to GB_Global struct
		/*const*/VC_Global* VC;		// pointer to VC_Global struct
		/*const*/chromosome* chroms;	// pointer to chromosomes' array
		/*const*/genlim* gene_lim;		// pointer to gene_lim genlim array (useful for bondaries defined for each gene)
		atom* atoms;				// pointer to atoms' array
		resid* residue;				// pointer to residues' array
		/*const*/gridpoint* cleftgrid;	// pointer to gridpoints' array (defining the total search space of the simulation)
		std::vector<float> 			Vectorized_Chromosome(chromosome* chrom);
		std::vector<float>			Vectorized_Cartesian_Coordinates(int chrom_index);
};

/*****************************************\
			RandomProjections
\*****************************************/
class RandomProjectedNeighborsAndDensities
{
	friend class FastOPTICS;
	
	public:
		explicit RandomProjectedNeighborsAndDensities(std::vector< std::pair< chromosome*,std::vector<float> > >&, int, FastOPTICS*); // Constructor (publicly called from FlexAID *_cluster.cxx)
	
	private:
		// private attributes
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
		void								getInverseDensities(std::vector<float> &);
		void								getNeighbors(std::vector< std::vector< int > > &);
		std::vector<float> 					Randomized_InternalCoord_Vector();
        std::vector<float>                  Randomized_CartesianCoord_Vector();
        std::vector<float>                  Randomly_Selected_Chromosome();
        bool 								accept_intraset_distance(std::vector<int>);
        void								output_projected_distance(char* end_strfile, char* tmp_end_strfile);
};
#endif
