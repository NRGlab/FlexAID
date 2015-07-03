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
		BindingMode();

		void 	add_Pose(chromosome*);
		double 	compute_energy();
		double 	compute_entropy();
		double 	compute_enthalpy();
		chromosome* elect_representative();

 	
 	protected:
		std::vector<Pose> Poses;
		BindingPopulation* Population; // used to access the BindingPopulation
		// void update_representative(); 	// force-update the BindingMode representative
		// chromosome* representative;
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
		void add_BindingMode(BindingMode&);
	private:
		std::vector< BindingMode > BindingModes;
};
#endif