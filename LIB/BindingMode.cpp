#include "BindingMode.h"

/*****************************************\
			BindingPopulation  
\*****************************************/

BindingPopulation::BindingPopulation(unsigned int temp) : Temperature(temp)
{
	this->PartitionFunction = 0.0;
}

void BindingPopulation::add_BindingMode(BindingMode mode)
{
	for(std::vector<Pose>::iterator pose = mode.Poses.begin(); pose != mode.Poses.end(); ++pose)
	{
		this->PartitionFunction += pose->boltzmann_weight;
	}

	this->BindingModes.push_back(mode);
	this->Entropize();
}

void BindingPopulation::Entropize() { std::sort(this->BindingModes.begin(), this->BindingModes.end(), BindingPopulation::EnergyComparator::EnergyComparator()); }

int BindingPopulation::get_BindingModes_size() { return this->BindingModes.size(); }
/*****************************************\
			  BindingMode
\*****************************************/
// public constructor
BindingMode::BindingMode(BindingPopulation* pop) : Population(pop) 
{}

// public method for pose addition
void BindingMode::add_Pose(chromosome* chrom, int chrom_index)
{
	Pose::Pose pose(chrom, chrom_index, this->Population->Temperature);
	this->Poses.push_back(pose);
}

double BindingMode::compute_enthalpy() const
{
	double enthalpy = 0.0;
	// compute enthalpy
	for(std::vector<Pose>::const_iterator pose = this->Poses.begin(); pose != this->Poses.end(); ++pose)
	{
		double boltzmann_prob = pose->boltzmann_weight / this->Population->PartitionFunction;
		enthalpy += boltzmann_prob * pose->CF;
	}
	return enthalpy;
}

double BindingMode::compute_entropy() const
{ 
	double entropy = 0.0;
	// compute entropy
	for(std::vector<Pose>::const_iterator pose = this->Poses.begin(); pose != this->Poses.end(); ++pose)
	{
		double boltzmann_prob = pose->boltzmann_weight / this->Population->PartitionFunction;
		entropy += boltzmann_prob * log(boltzmann_prob);
	}
	return entropy;
}

double BindingMode::compute_energy() const
{ 
	return ( this->compute_enthalpy() + ( this->Population->Temperature * this->compute_entropy() ) );
}


/*****************************************\
				  Pose
\*****************************************/

Pose::Pose(chromosome* chrom, int index, uint temperature) : chrom(chrom), chrom_index(index), CF(chrom->app_evalue)
{
	this->boltzmann_weight = pow( E, ((-1.0) * (1/static_cast<double>(temperature)) * chrom->app_evalue) );
}