#ifndef BINDINGMODE_H
#define BINDINGMODE_H
#include "gaboom.h"
/*****************************************\
				  Pose
\*****************************************/
class Pose
{
	friend class BindingMode;

	public:

	protected:
		chromosome* chrom;
		double boltzmann_prob;
};

/*****************************************\
			  BindingMode
\*****************************************/
class BindingMode // aggregation of poses (Cluster)
{
	friend class Clusters;

	public:
		void 	add_pose();
		double 	get_energy();
		double 	get_entropy();
		double 	get_enthalpy();

	protected:
		void update_representative();
		void update_enthalpy();
		void update_entropy();
		void update_energy(uint temperature);

		chromosome* representative;
		std::vector<Pose> poses;
		double enthalpy;
		double entropy;
		double energy;
};
// * add to -> BindingMode.cpp
// public 
void BindingMode::add_pose(chromosome* pose)
{

}
double get_enthalpy() { return this->enthalpy; }
double get_entropy()  { return this->entropy; }
double get_energy()   { return this->energy; }
//private 

/*****************************************\
				Clusters  
\*****************************************/
class Clusters
{

};
#endif