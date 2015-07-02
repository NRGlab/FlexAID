#include "BindingMode.h"
/*****************************************\
			BindingPopulation  
\*****************************************/
BindingPopulation::BindingPopulation(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chroms, genelim* gene_lim, atom* atoms, resid* resid, gridpoint* grid, int num_chrom) : FA(FA), GB(GB), VC(VC), chromosomes(chroms), gene_lim(gene_lim), atoms(atoms), residue(resid), cleftgrid(grid), Temperature(FA->temperature), nChrom(num_chrom)
{
	this->PartitionFunction = 0.0;
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
double get_enthalpy() { return this->enthalpy; }
double get_entropy()  { return this->entropy; }
double get_energy()   { return this->energy; }
//private 
/*****************************************\
				  Pose
\*****************************************/
Pose::Pose(chromosome* chrom, uint temperature) : chrom(chrom)
{
	this->boltzmann_weight = pow( E, ((-1.0) * (1/static_cast<double>(temperature)) * this->chrom->app_evalue) );
}