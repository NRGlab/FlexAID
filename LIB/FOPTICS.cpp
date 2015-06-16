#include "FOPTICS.h"

/*****************************************\
				FastOPTICS
\*****************************************/
// Constructor
FastOPTICS::FastOPTICS(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom)
{
	FastOPTICS::Initialize(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,num_chrom);
	
	// vector of point indexes 
	std::vector< int > ptInd;
	for(int k = 0; k < this->N; ++k) ptInd.push_back(k);
	
	RandomProjectedNeighborsAndDensities::RandomProjectedNeighborsAndDensities MultiPartition(this->points, this->minPts);
	MultiPartition.computeSetBounds(ptInd);

}

FastOPTICS::Initialize(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chroms, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom)
{
	// FlexAID
	int this->N = nChrom;
	FA_Global* FastOPTICS::FA = FA;
	GB_Global* FastOPTICS::GB = GB;
	VC_Global* FastOPTICS::VC = VC;
	chromosome* FastOPTICS::chroms = chroms;
	int FastOPTICS::minPts = (int)(0.01*this->N - 0.5);
	int FastOPTICS::nDimensions = this->FA->npar + 2; // 3 Dim for first gene (translational) + 1 Dim per gene = nGenes + 2
	gene_lim* FastOPTICS::gene_lim = gene_lim;
	atom* FastOPTICS::atoms = atoms;
	resid* FastOPTICS::residue = residue;
	gridpoint* FastOPTICS::cleftgrid = cleftgrid;
	// initialization from chromosome* chrom (into a pointer or into vector<>)

	// FastOPTICS
	int FastOPTICS::iOrder = 0;
	std::vector< int > 		this->FastOPTICS::order(this->N, 0);
	std::vector< float > 	this->FastOPTICS::reachDist(this->N, UNDEFINED_DIST);
	std::vector< bool > 	this->FastOPTICS::processed(this->N, false);
	std::vector< float > 	this->FastOPTICS::inverseDensities(this->N, 0.0f);
	this->order.reserve(this->N);
	this->reachDist.reserve(this->N);
	this->processed.reserve(this->N);
	this->inverseDensities.reserve(this->N);

	std::vector< std::pair< chromosome*,std::vector<float> > > this->FastOPTICS::points;
	this->points.reserve(this->N);
	std::vector< std::vector< chromosome* > > this->FastOPTICS::neighs(this->N);
	this->neighs.reserve(this->N);
	
	for(int i = 0; i < this->N; ++i)
	{
		// need to transform the chromosomes into vector f
		std::vector<float> vChrom(FastOPTICS::vectorize_chrom(&chroms[i])); // Copy constructor
		if(vChrom.size() == this->nDimensions)
			// std::pair<first, second> is pushed to this->points[]
			// 	first  -> chromosome* pChrom (pointer to chromosome)
			// 	second -> vChrom[FA->npar+n] (vectorized chromosome)
			this->points.push_back( std::make_pair( (chromosome*)&chroms[i], vChrom) ) ;
	}
}

std::vector<float> FastOPTICS::vectorize_chromosome(chromsome* chrom)
{
	std::vector<float> vChrom(this->nDimensions, 0.0f);
	// getting 
	for(int j = 0; j < this->nDimensions, ++j)
	{
		if(j == 0) //  building the first 3 comp. from genes[0] which are CartCoord x,y,z
			for(int i = 0; i < 3; ++i)
				vChrom[i] = (float) this->cleftgrid[(unsigned int)(*chrom).genes[j].to_ic].coor[i] - this->FA->ori[i];
		else
		{
			// j+2 is used from {j = 1 to N} to build further comp. of genes[j]
			vChrom[j+2] = (float) genetoic(&gene_lim[j], (*chrom).genes[j].to_int32);
		}
	}
	return vChrom;
}

/*****************************************\
			RandomProjections
\*****************************************/
// Constructor
RandomProjectedNeighborsAndDensities::RandomProjectedNeighborsAndDensities(std::vector< std::pair< chromosome*,std::vector<float> > >& points, int minSplitSize)
{
	this->minSplitSize = minSplitSize;

	if( !points || points.empty() )
	{
		this->N = 0;
		this->nDimensions = 0;
		this->nProject1D = 0;
		return;
	}
	else
	{
		this->N = this->points.size();
		this->nDimensions = (this->points[0].second).size();
		this->nPointsSetSplits = (int) this->logOProjectionConstant * log(this->N * this->nDimensions + 1)/log(2);
		this->nProject1D = (int) this->logOProjectionConstant * log(this->N * this->nDimensions + 1)/log(2);
		// line below calls copy-constructor
		std::vector< std::pair< chromosome*,std::vector<float> > > this->RandomProjectedNeighborsAndDensities::points(points);
		// construct splitsets && projectPoints by calling vector constructor
		// default constructor for empty vector for splitsets
		std::vector< std::vector< int > > this->RandomProjectedNeighborsAndDensities::splitsets;
		// default constructor for size this->nProject1D for projectedPoints[nProject1D][N]
		std::vector< std::vector<float> > this->RandomProjectedNeighborsAndDensities::projectedPoints;
		// std::vector< std::vector<float> > RandomProjectedNeighborsAndDensities::projectedPoints(this->nProject1D);
	}
}

void RandomProjectedNeighborsAndDensities::computeSetBounds(std::vector< int > & ptList)
{
	std::vector< std::vector<float> > tempProj;
	// projecting points
	for(int j = 0; j<this->nProject1D; ++j)
	{
		std::vector<float> currentRp(this->nDimensions,0.0f);
		for(int i=0; i<this->nDimensions; ++i)
		{
			if(i==0) // First 3Dim are defined by the grid on genes[0]
			{

			}
			else // other Dim are defined by subsequent genes[i] ? (i > 0)
			{

			}
		}
	}
}