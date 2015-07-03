#ifndef BINDINGMODE_H
#define BINDINGMODE_H

#include "gaboom.h"
#include "boinc.h"

/*****************************************\
				  Pose
\*****************************************/
struct Pose
{
	// public constructor :
	Pose(chromosome* chrom, uint temperature);
	// public (default behavior when struct is used instead of class)
	chromosome* chrom;
	double CF;
	double boltzmann_weight;
};
/*****************************************\
			  BindingMode
\*****************************************/
class BindingPopulation; // forward-declaration in order to access BindingPopulation* Population pointer
class BindingMode // aggregation of poses (Cluster)
{
	friend class BindingPopulation;
	
	public:
		BindingMode(BindingPopulation*);

		void 	add_Pose(chromosome*);
		double 	compute_energy() const;
		double 	compute_entropy() const;
		double 	compute_enthalpy() const;
		chromosome* elect_representative();

 	protected:
		std::vector<Pose> Poses;
		BindingPopulation* Population; // used to access the BindingPopulation
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
		// public constructor to be called once
		BindingPopulation(unsigned int); 
		void add_BindingMode(BindingMode);
	
	private:
		std::vector< BindingMode > BindingModes;
		
		void Entropize(); // Sort BindinModes according to their observation frequency
		
		struct EnergyComparator
		{
			inline bool operator() ( const BindingMode& BindingMode1, const BindingMode& BindingMode2 )
			{
				double energy1 = BindingMode1.compute_energy();
				double energy2 = BindingMode2.compute_energy();
				return (energy1 < energy2);
			}
		};
};
#endif