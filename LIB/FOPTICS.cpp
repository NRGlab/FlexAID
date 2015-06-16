#include "FOPTICS.h"

/*****************************************\
				FastOPTICS
\*****************************************/
// Constructor
FastOPTICS::FastOPTICS(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom)
{
	this->FastOPTICS::Initialize(FA,GB,VC,chrom,gene_lim,atoms,residue,cleftgrid,num_chrom);
	
	// vector of point indexes 
	std::vector< int > ptInd;
	for(int k = 0; k < this->N; ++k) ptInd.push_back(k);
	
	RandomProjectedNeighborsAndDensities::RandomProjectedNeighborsAndDensities MultiPartition(this->points, this->minPoints, this);
	MultiPartition.computeSetBounds(ptInd);

}

void FastOPTICS::Initialize(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chroms, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int num_chrom)
{
	// FlexAID
	this->N = num_chrom;
	this->minPoints = (int)(0.01*this->N - 0.5);
	FastOPTICS::FA = FA;
	FastOPTICS::GB = GB;
	FastOPTICS::VC = VC;
	this->nDimensions = this->FA->npar + 2; // 3 Dim for first gene (translational) + 1 Dim per gene = nGenes + 2
	FastOPTICS::chroms = chroms;
	FastOPTICS::gene_lim = gene_lim;
	FastOPTICS::atoms = atoms;
	FastOPTICS::residue = residue;
	FastOPTICS::cleftgrid = cleftgrid;
	// initialization from chromosome* chrom (into a pointer or into vector<>)

	// FastOPTICS
	FastOPTICS::iOrder = 0;
	// std::vector< int > 		this->FastOPTICS::order(this->N, 0);
	// std::vector< float > 	this->FastOPTICS::reachDist(this->N, UNDEFINED_DIST);
	// std::vector< bool > 	this->FastOPTICS::processed(this->N, false);
	// std::vector< float > 	this->FastOPTICS::inverseDensities(this->N, 0.0f);
	this->order.reserve(this->N);
	this->reachDist.reserve(this->N);
	this->processed.reserve(this->N);
	this->inverseDensities.reserve(this->N);

	// std::vector< std::pair< chromosome*,std::vector<float> > > this->FastOPTICS::points;
	this->points.reserve(this->N);
	// std::vector< std::vector< chromosome* > > this->FastOPTICS::neighs(this->N);
	this->neighs.reserve(this->N);
	
	for(int i = 0; i < this->N; ++i)
	{
		// need to transform the chromosomes into vector f
		std::vector<float> vChrom(this->FastOPTICS::Vectorized_Chromosome(&chroms[i])); // Copy constructor
		if(vChrom.size() == this->nDimensions)
			// std::pair<first, second> is pushed to this->points[]
			// 	first  -> chromosome* pChrom (pointer to chromosome)
			// 	second -> vChrom[FA->npar+n] (vectorized chromosome)
			this->points.push_back( std::make_pair( (chromosome*)&chroms[i], vChrom) ) ;
	}
}

std::vector<float> FastOPTICS::Vectorized_Chromosome(chromosome* chrom)
{
	std::vector<float> vChrom(this->nDimensions, 0.0f);
	// getting 
	for(int j = 0; j < this->nDimensions; ++j)
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
RandomProjectedNeighborsAndDensities::RandomProjectedNeighborsAndDensities(std::vector< std::pair< chromosome*,std::vector<float> > >& points, int minSplitSize, FastOPTICS* top)
{
	this->top = top;
	this->minSplitSize = minSplitSize;

	if( points.empty() )
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
		this->points = points;
	}
}

void RandomProjectedNeighborsAndDensities::computeSetBounds(std::vector< int > & ptList)
{
	std::vector< std::vector<float> > tempProj;
	// perform projection of points
	for(int j = 0; j<this->nProject1D; ++j)
	{
		std::vector<float> currentRp = this->Randomized_Normalized_Vector();

		int k = 0;
		std::vector<int>::iterator it = ptList.begin();
		while(it != ptList.end())
		{
			float sum = 0.0f;
			std::vector<float>::iterator vecPt = this->points[(*it)++].second.begin();
			std::vector<float>::iterator currPro = this->projectedPoints[j].begin();
			for(int m = 0; m < this->nDimensions; ++m)
			{
				sum += currentRp[m] * vecPt[m];
			}
			currPro[k] = sum;
			++k;
		}
	}

	// Split Points Set
	std::vector<int> projInd;
	projInd.reserve(this->nProject1D);
	for(int j = 0; j < this->nProject1D; ++j) 
		projInd.push_back(j);
	for(int avgP = 0; avgP < this->nPointsSetSplits; avgP++)
	{
		// Shuffle projections
		for(int i = 0; i < this->nProject1D; ++i)
			tempProj.push_back(this->projectedPoints[i]);
		std::random_shuffle(projInd.begin(), projInd.end());
		std::vector<int>::iterator it = projInd.begin();
		int i = 0;
		while(it != projInd.end())
		{
			int cind = (*it)++;
			this->projectedPoints[cind] = tempProj[i];
			i++;
		}

		//split points set
		int nPoints = ptList.size();
		std::vector<int> ind(nPoints);
		ind.reserve(nPoints);
		for(int l = 0; l < nPoints; ++l)
			ind[l] = l;
		this->SplitUpNoSort(ind,0);
	}
}

void RandomProjectedNeighborsAndDensities::SplitUpNoSort(std::vector< int >& ind, int dim)
{
	int nElements = ind.size();
	dim = dim % this->nProject1D;
	std::vector<float>::iterator tProj = this->projectedPoints[dim].begin();
	int splitPos;


}

std::vector<float> RandomProjectedNeighborsAndDensities::Randomized_Normalized_Vector()
{
	std::vector<float> vChrom(this->nDimensions, 0.0f);
	
	float sum = 0.0f;
	for(int j = 0; j < this->nDimensions; ++j)
	{
		double intGene = this->Dice();
		double doubleGene = genetoic(&this->top->gene_lim[j],intGene);
		if(j == 0) //  building the first 3 comp. from genes[0] which are CartCoord x,y,z
			for(int i = 0; i < 3; ++i)
			{
				vChrom[i] = (float) this->top->cleftgrid[(unsigned int)doubleGene].coor[i] - this->top->FA->ori[i];
				sum += vChrom[i]*vChrom[i];
			}
		else
		{
			// j+2 is used from {j = 1 to N} to build further comp. of genes[j]
			vChrom[j+2] = (float) doubleGene;
			sum += vChrom[j+2]*vChrom[j+2];
		}
	}
	sum = sqrt(sum);
	for(int k = 0; k < this->nDimensions; ++k) vChrom[k]/=sum;

	return vChrom;
}

int RandomProjectedNeighborsAndDensities::Dice()
{
	unsigned int tt = static_cast<unsigned int>(time(0));
	srand(tt);
	RNGType rng(tt);
	boost::uniform_int<> one_to_max_int32( 0, MAX_RANDOM_VALUE );
	boost::variate_generator< RNGType, boost::uniform_int<> > dice(rng, one_to_max_int32);
	return dice();
}