#ifndef BINDINGMODE_H
#define BINDINGMODE_H

#include "gaboom.h"
#include "boinc.h"
#include <utility>
#include <queue>

/*****************************************\
				  Pose
\*****************************************/
struct Pose
{
	// public constructor :
	Pose(chromosome* chrom, uint temperature);
	// public (default behavior when struct is used instead of class)
	chromosome* chrom;
	double boltzmann_weight;
};
/*****************************************\
			  BindingMode
\*****************************************/
class BindingPopulation; // forward-declaration in order to access BindingPopulation* Population pointer
class BindingMode // aggregation of poses (Cluster)
{
	public:
		BindingMode();

		void 	add_Pose(chromosome*);
		double 	get_energy();
		double 	get_entropy();
		double 	get_enthalpy();

	protected:
		void update_representative(); 	// force-update the BindingMode representative
		void update_enthalpy();			// recalculate BindingMode enthalpy
		void update_entropy();			// recalculate BingingMode entropy
		void update_energy();

		chromosome* representative;
		std::vector<Pose> Poses;
		double enthalpy;
		double entropy;
		double energy;

	private:
		BindingPopulation* Population; // used to access the BindingPopulation
};

/*****************************************\
			BindingPopulation  
\*****************************************/
class BindingPopulation
{
	public:
		// public constructor to be called once
		BindingPopulation(FA_Global*, GB_Global*, VC_Global*, chromosome*, genlim*, atom*, resid*, gridpoint*, int); 
		void add_BindingMode(BindingMode);
	private:
		std::vector< BindingMode > BindingModes;

	protected:
		// Temperature is used for energy calculations of BindingModes
		int nChrom;
		unsigned int Temperature;
		double PartitionFunction;
		// FlexAID data structure references' pointers
		FA_Global* FA;
		GB_Global* GB;
		VC_Global* VC;
		atom* atoms;
		chromosome* chromosomes;
		genlim* gene_lim;
		resid* residue;
		gridpoint* cleftgrid;
};
#endif