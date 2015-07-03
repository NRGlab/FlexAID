#include "BindingMode.h"
/*****************************************\
			BindingPopulation  
\*****************************************/
BindingPopulation::BindingPopulation(unsigned int temp) : Temperature(temp) 
{
	this->PartitionFunction = 0.0;
}

void BindingPopulation::add_BindingMode(BindingMode& mode)
{
	for(std::vector<Pose>::iterator pose = mode.Poses.begin(); pose != mode.Poses.end(); ++pose)
		this->PartitionFunction += pose->boltzmann_weight;

	this->BindingModes.push_back(mode);
}
/*****************************************\
			  BindingMode
\*****************************************/
// public constructor
BindingMode::BindingMode() 
{}
// public method for pose addition
void BindingMode::add_Pose(chromosome* chrom)
{
	Pose::Pose pose(chrom, this->Population->Temperature);
	this->Poses.push_back(pose);
}
double BindingMode::compute_enthalpy()
{
	double enthalpy = 0.0;
	// compute enthalpy
	for(std::vector<Pose>::iterator pose = this->Poses.begin(); pose != this->Poses.end(); ++pose)
	{
		double boltzmann_prob = pose->boltzmann_weight / this->Population->PartitionFunction;
		enthalpy += boltzmann_prob * pose->CF;
	}
	return enthalpy;
}
double BindingMode::compute_entropy()
{ 
	double entropy = 0.0;
	// compute entropy
	for(std::vector<Pose>::iterator pose = this->Poses.begin(); pose != this->Poses.end(); ++pose)
	{
		double boltzmann_prob = pose->boltzmann_weight / this->Population->PartitionFunction;
		entropy += boltzmann_prob * log(boltzmann_prob);
	}
	return entropy;
}
double BindingMode::compute_energy()
{ 
	return this->compute_enthalpy() + this->Population->Temperature * this->compute_entropy();
}
//private 
/*****************************************\
				  Pose
\*****************************************/
Pose::Pose(chromosome* chrom, uint temperature) : chrom(chrom), CF(chrom->app_evalue)
{
	this->boltzmann_weight = pow( E, ((-1.0) * (1/static_cast<double>(temperature)) * chrom->app_evalue) );
}