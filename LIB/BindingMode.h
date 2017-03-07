#ifndef BINDINGMODE_H
#define BINDINGMODE_H

#include "gaboom.h"
#include "boinc.h"
#include <ctime>

//#define UNDEFINED_DIST FLT_MAX // Defined in FOPTICS as > than +INF
#define UNDEFINED_DIST -0.1f // Defined in FOPTICS as > than +INF
#define isUndefinedDist(a) ((a - UNDEFINED_DIST) <= FLT_EPSILON)

int roll_die();

class BindingPopulation; // forward-declaration in order to access BindingPopulation* Population pointer
/*****************************************\
				  Pose
\*****************************************/
struct Pose
{
	// public constructor :
	Pose(chromosome* chrom, int chrom_index, int order, float dist, uint temperature, std::vector<float>);
	~Pose();
	// public (default behavior when struct is used instead of class)
	bool processed;
	int chrom_index;
	int order;
	float reachDist;
	chromosome* chrom;
	double CF;
	double CFdS;
	double boltzmann_weight;
	std::vector<float> vPose;
	inline bool const operator< (const Pose& rhs);
	inline bool const operator> (const Pose& rhs);
	inline bool const operator==(const Pose& rhs);
};
// publicly accessible 2 argumentss operator== (lhs, rhs)
inline bool const operator==(const Pose& lhs, const Pose& rhs);

struct PoseClassifier
{
   inline bool operator() ( const Pose& Pose1, const Pose& Pose2 )
   {
       	if(Pose1.order < Pose2.order) return true;
       	else if(Pose1.order > Pose2.order) return false;
		
		if(Pose1.reachDist < Pose2.reachDist) return true;
		else if(Pose1.reachDist > Pose2.reachDist) return false;
		
		if(Pose1.chrom_index < Pose2.chrom_index) return true;
		else if(Pose1.chrom_index > Pose2.chrom_index) return false;
		
		return false;
   }
};

struct PoseRanker
{
	inline bool operator() ( const Pose& Pose1, const Pose& Pose2 )
   {
   		if(Pose1.CFdS < Pose2.CFdS) return true;
   		else if(Pose1.CFdS > Pose2.CFdS) return false; 

       	if(Pose1.CF < Pose2.CF) return true;
       	else if(Pose1.CF > Pose2.CF) return false;
		
		if(Pose1.boltzmann_weight > Pose2.boltzmann_weight) return true;
		else if(Pose1.boltzmann_weight < Pose2.boltzmann_weight) return false;
		
		if(Pose1.chrom_index < Pose2.chrom_index) return true;
		else if(Pose1.chrom_index > Pose2.chrom_index) return false;
		
		return false;
   }
};
/*****************************************\
			  BindingMode
\*****************************************/
class BindingMode // aggregation of poses (Cluster)
{
	friend class BindingPopulation;
	friend class FastOPTICS;
	friend class ColonyEnergy;

	public:
												// public constructor (explicitely needs a pointer to a BindingPopulation of type BindingPopulation*)
		explicit 								BindingMode(BindingPopulation*);

			void								add_Pose(Pose&);
												//	tells if a Pose can be added to a BidingMode
			bool 								isPoseAggregable(const Pose& pose) const;
			bool								isPoseInBindingMode(int chrom_index) const;
												//	tells if the BindingMode contains one *homogenic* population
			bool 								isHomogenic() const;
			void                                clear_Poses();
			int									get_BindingMode_size() const;
			float 								compute_distance(const Pose& pose1, const Pose& pose2) const;
			double								compute_energy() const;
			double								compute_entropy() const;
			double								compute_enthalpy() const;
			std::vector<Pose>::const_iterator 	elect_Representative(bool useCentroid) const;
			inline bool const 					operator< (const BindingMode& rhs);
			inline bool const 					operator==(const BindingMode& rhs);


 	protected:
		std::vector<Pose> Poses;
		BindingPopulation* Population; // used to access the BindingPopulation

		bool 	isValid;
		void	set_energy();

	private:
		// private attributes
		double 	energy;
		// private methods 
		std::vector<float>	compute_centroid() const;	// return the centroid as a N dimensional vector of cartesian coordinates
		double 				compute_partition_function() const; // used in compute_centroid
		void 				output_BindingMode(int num_result, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, int minPoints);
		void				output_dynamic_BindingMode(int nBindingMode, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, int minPoints);
};
// publicly accessible 2 argumentss operator== (lhs, rhs)
inline bool const operator==(const BindingMode& lhs, const BindingMode& rhs);
/*****************************************\
			BindingPopulation  
\*****************************************/
class BindingPopulation
{
	friend class BindingMode;
    friend class FastOPTICS;
    friend class ColonyEnergy;
    
	public:
		// Temperature is used for energy calculations of BindingModes
		unsigned int Temperature;
		// Explicit constructor to be called (it will handle the PartitionFunction)
		explicit 	BindingPopulation(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom);
			// 	add new binding mode to population
		void				add_BindingMode(BindingMode&);
			// Classify_BindingModes merges similar BindingModes together
		void 				Classify_BindingModes();
			// used to compute the distance between 2 poses (same as BindingMode::compute_distance)
		float 				compute_distance(const Pose& pose1, const Pose& pose2) const;
		float 				compute_distance(std::vector<Pose>::const_iterator,std::vector<Pose>::const_iterator) const;
		float 				compute_vec_distance(std::vector<float>, std::vector<float>) const;
			// 	return the number of BindinMonde (size getter)
		int					get_Population_size();
			// 	output BindingMode up to nResults results
		void				output_Population(int nResults, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp, int minPoints);
		std::vector<float> 	Vectorized_Chromosome(chromosome* chrom);
		std::vector<float>	Vectorized_Cartesian_Coordinates(int chrom_index);

		std::vector< Pose > Poses; 					// Poses container

	protected:
		double PartitionFunction;	// sum of all Boltzmann_weight
		int nChroms;				// n_chrom_snapshot input to clustergin function
		int nDimensions;

		// FlexAID pointers
		FA_Global* 	FA;			// pointer to FA_Global struct
		GB_Global* 	GB;			// pointer to GB_Global struct
		VC_Global* 	VC;			// pointer to VC_Global struct
		chromosome* chroms;		// pointer to chromosomes' array
		genlim* gene_lim;		// pointer to gene_lim genlim array (useful for bondaries defined for each gene)
		atom* atoms;			// pointer to atoms' array
		resid* residue;			// pointer to residues' array
		gridpoint* cleftgrid;	// pointer to gridpoints' array (defining the total search space of the simulation)

	        // Sort BindinModes according to their energy score (∆G = ∆H - T∆S)
        void 	Entropize();			// entropize BindingModes according to their current computed energy
	
	private:
		std::vector< BindingMode > 	BindingModes;	// BindingMode container
		
			// merges two existing BindingModes
		bool 	merge_BindingModes(std::vector< BindingMode >::iterator mode1, std::vector< BindingMode >::iterator mode2);
			// removes a BindingMode from the Ppopulation (THIS WILL NOT AFFACT THE POPULATION PARTITION FUNCTION)
		void 	remove_BindingMode(std::vector<BindingMode>::iterator mode);
			// removes all BindingModes where isValid == false (calls remove_BidningMode)
		void 	remove_invalid_BindingModes();
		struct EnergyComparator
		{
			inline bool operator() ( const BindingMode& BindingMode1, const BindingMode& BindingMode2 )
			{
				return (BindingMode1.compute_energy() < BindingMode2.compute_energy());
			}
		};
};
#endif
