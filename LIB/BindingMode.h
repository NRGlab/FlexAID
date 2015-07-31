#ifndef BINDINGMODE_H
#define BINDINGMODE_H

#include "gaboom.h"
#include "boinc.h"

// Float comparators
// bool definitelyGreaterThan(float a, float b, float epsilon);
// bool definitelyLessThan(float a, float b, float epsilon);

// bool definitelyGreaterThan(float a, float b, float epsilon)
// {
//     return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }

// bool definitelyLessThan(float a, float b, float epsilon)
// {
//     return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }
/*****************************************\
				  Pose
\*****************************************/
struct Pose
{
	// public constructor :
	Pose(chromosome* chrom, int chrom_index, int order, float dist, uint temperature);
	// public (default behavior when struct is used instead of class)
	int chrom_index;
	int order;
	float reachDist;
	chromosome* chrom;
	double CF;
	double boltzmann_weight;
	inline bool const operator< (const Pose& rhs);
};

struct PoseClassifier
{
   inline bool operator() ( const Pose& Pose1, const Pose& Pose2 )
   {
       	if(Pose1.order < Pose2.order) return true;
       	else if(Pose1.order > Pose2.order) return false;
		else if(Pose1.reachDist < Pose2.reachDist) return true;
		else if(Pose1.reachDist > Pose2.reachDist) return false;
		else if(Pose1.CF < Pose2.CF) return true;
		else if(Pose1.CF > Pose2.CF) return false;
		else if(Pose1.chrom_index < Pose2.chrom_index) return true;
		else if(Pose1.chrom_index > Pose2.chrom_index) return false;
		else return false;
   }
};
/*****************************************\
			  BindingMode
\*****************************************/
class BindingPopulation; // forward-declaration in order to access BindingPopulation* Population pointer
class BindingMode // aggregation of poses (Cluster)
{
	friend class BindingPopulation;
	
	public:
		explicit 	BindingMode(BindingPopulation*); // public constructor (explicitely needs a pointer to a BindingPopulation of type BindingPopulation*)

			void		add_Pose(Pose&);
			void		clear_Poses();
			int			get_BindingMode_size() const;
			double		compute_energy() const;
			double		compute_entropy() const;
			double		compute_enthalpy() const;
			chromosome* elect_representative();
			inline bool const operator< (const BindingMode& rhs);

 	protected:
		std::vector<Pose> Poses;
		BindingPopulation* Population; // used to access the BindingPopulation

		void	set_energy();

	private:
		double energy;
};

/*****************************************\
			BindingPopulation  
\*****************************************/
class BindingPopulation
{
	public:
		// Temperature is used for energy calculations of BindingModes
		unsigned int Temperature;
		double PartitionFunction;
		
		explicit 	BindingPopulation(unsigned int);// public constructor (explicitely needs an int representative of Temperature)
			void	add_BindingMode(BindingMode&); 	// add new binding mode to population
			int		get_BindingModes_size();		// return the number of BindinMonde (size getter)
	
		void 		Entropize(); 					// Sort BindinModes according to their observation frequency

	private:
		std::vector< BindingMode > BindingModes;	// BindingMode container
		
		struct EnergyComparator
		{
			inline bool operator() ( const BindingMode& BindingMode1, const BindingMode& BindingMode2 )
			{
				return (BindingMode1.energy < BindingMode2.energy);
				// return (BindingMode1.compute_energy() < BindingMode2.compute_energy());
			}
		};
};
#endif