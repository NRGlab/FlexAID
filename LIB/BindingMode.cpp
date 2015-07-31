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
	{
		this->PartitionFunction += pose->boltzmann_weight;
	}
    mode.set_energy();
	this->BindingModes.push_back(mode);
	this->Entropize();
}

void BindingPopulation::Entropize()
{
	for(std::vector<BindingMode>::iterator it = this->BindingModes.begin(); it != this->BindingModes.end(); ++it)
	{
		it->set_energy();
	}
	std::sort(this->BindingModes.begin(), this->BindingModes.end(), BindingPopulation::EnergyComparator::EnergyComparator());
}

int BindingPopulation::get_BindingModes_size() { return this->BindingModes.size(); }
/*****************************************\
			  BindingMode
\*****************************************/
// public constructor
BindingMode::BindingMode(BindingPopulation* pop) : Population(pop), energy(0.0)
{
}

// public method for pose addition
void BindingMode::add_Pose(Pose& pose)
{
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

int BindingMode::get_BindingMode_size() const { return this->Poses.size(); }

void BindingMode::clear_Poses() { this->Poses.clear(); }

void BindingMode::set_energy()
{
	this->energy = this->compute_energy();
}
inline bool const BindingMode::operator< (const BindingMode& rhs) { return this->compute_energy() < rhs.compute_energy(); }
/*****************************************\
				  Pose
\*****************************************/
Pose::Pose(chromosome* chrom, int index, int iorder, float dist, uint temperature) : chrom(chrom), order(iorder), chrom_index(index), reachDist(dist), CF(chrom->app_evalue)
{
	this->boltzmann_weight = pow( E, ((-1.0) * (1/static_cast<double>(temperature)) * chrom->app_evalue) );
}
inline bool const Pose::operator< (const Pose& rhs)
{
	if(this->order < rhs.order) return true;
   	else if(this->order > rhs.order) return false;
	else if(this->reachDist < rhs.reachDist) return true;
	else if(this->reachDist > rhs.reachDist) return false;
	else if(this->CF < rhs.CF) return true;
	else if(this->CF > rhs.CF) return false;
	else if(this->chrom_index < rhs.chrom_index) return true;
	else if(this->chrom_index > rhs.chrom_index) return false;
	else return false;
}