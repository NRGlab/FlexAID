#include "FOPTICS.h"

/*****************************************\
				FastOPTICS
\*****************************************/
// Constructor
FastOPTICS::FastOPTICS(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom)
{
	this->FastOPTICS::Initialize(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,nChrom);

	RandomProjectedNeighborsAndDensities::RandomProjectedNeighborsAndDensities algo(, minPts);


}

FastOPTICS::Initialize(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom)
{
	// FlexAID
	this->N = nChrom;
	this->FA = FA;
	this->GB = GB;
	this->VC = VC;
	// initialization from chromosome* chrom (into a pointer or into vector<>)
	this->gene_lim = gene_lim;
	this->atoms = atoms;
	this->residue = residue;
	this->cleftgrid = cleftgrid;

	// FastOPTICS
	this->iOrder = 0;

	std::vector< int > FastOPTICS::order(this->N);
	std::vector< double > FastOPTICS::reachDist(this->N);
	std::vector< bool > FastOPTICS::processed(this->N);
	std::vector< double > FastOPTICS::inverseDensities(this->N);
	std::vector< chromosome* > FastOPTICS::chroms(this->N);
	std::vector< std::vector< chromosome* > > FastOPTICS::neighs(this->N);
	
	for(int i = 0; i < this->N; ++i)
	{
		this->order.push_back(0);
		this->reachDist.push_back(UNDEFINED_DIST);
		this->processed.push_back(false);
		this->inverseDensities.push_back(0.0);
		this->chroms.push_back((chromosome*)&chrom[i]);
	}
}