#include "FOPTICS.h"

ClusterOrdering::ClusterOrdering(int id, int predID, float reach) : objectID(id), predecessorID(predID), reachability(reach)
{
	// this->objectID = id;
	// this->predecessorID = predID;
	// this->reachability = reach;
}

inline bool const ClusterOrdering::operator==(const ClusterOrdering& rhs)
{
	if(*this == rhs)
		return true;
	if(typeid(rhs) != typeid(ClusterOrdering))
		return false;

	return ( this->objectID == rhs.objectID );
}

inline bool const ClusterOrdering::operator< (const ClusterOrdering& rhs)
{
	if(this->reachability > rhs.reachability)
			return false;
		else if(this->reachability < rhs.reachability)
			return true;
		if(this->objectID > rhs.objectID)
			return false;
		else if(this->objectID < rhs.objectID)
			return true;
		// if nothing else is true, return 0
		return false;
}
/*****************************************\
				FastOPTICS
\*****************************************/
// Constructor and Algorithm main+only call
int FastOPTICS::iOrder;
FastOPTICS::FastOPTICS(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, gridpoint* cleftgrid, int num_chrom, BindingPopulation& Population)
{	
// Declarations
	///////////////////////////////////////////////////////
	// Entropy
	this->Population = &Population;
	// FlexAID
	this->N = num_chrom;
//	this->minPoints = 5;
	this->minPoints = static_cast<int>( floor(this->N * 0.01) );
	this->FA = FA;
	this->GB = GB;
	this->VC = VC;
	this->cleftgrid = cleftgrid;
	this->nDimensions = this->FA->npar + 2; // 3 Dim for first gene (translational) + 1 Dim per gene = nGenes + 2
	this->chroms = chrom;
	this->gene_lim = gene_lim;

	// FastOPTICS
	FastOPTICS::iOrder = 0;
	this->order.reserve(this->N);
	this->reachDist.reserve(this->N);
	this->processed.reserve(this->N);
	this->inverseDensities.reserve(this->N);
	this->points.reserve(this->N);
	this->neighbors.reserve(this->N);
	
	for(int i = 0; i < this->N; ++i)
	{
		// need to transform the chromosomes into vector f
		std::vector<float> vChrom(this->FastOPTICS::Vectorized_Chromosome(&chrom[i])); // Copy constructor
		if(vChrom.size() == this->nDimensions)
		{
			// std::pair<first, second> is pushed to this->points[]
			// 	first  -> chromosome* pChrom (pointer to chromosome)
			// 	second -> vChrom[FA->npar+n] (vectorized chromosome)
			this->order.push_back(-1);
			this->reachDist.push_back(FLT_MAX);
			this->processed.push_back(false);
			this->inverseDensities.push_back(0.0f);
			this->points.push_back( std::make_pair( (chromosome*)&chrom[i], vChrom) ) ;
		}
	}

}

void FastOPTICS::Execute_FastOPTICS()
{
	// vector of point indexes 
	std::vector< int > ptInd;
    ptInd.reserve(this->N);
    for(int k = 0; k < this->N; ++k)
    {
        ptInd.push_back(k);
    }
	
	// Build object, compute projections, density estimates and density neighborhoods (in serial order of function calls below)
	RandomProjectedNeighborsAndDensities::RandomProjectedNeighborsAndDensities MultiPartition(this->points, this->minPoints, this);
	MultiPartition.RandomProjectedNeighborsAndDensities::computeSetBounds(ptInd);
	MultiPartition.RandomProjectedNeighborsAndDensities::getInverseDensities(this->inverseDensities);
	MultiPartition.RandomProjectedNeighborsAndDensities::getNeighbors(this->neighbors);

	// Compute OPTICS ordering
	for(int ipt = 0; ipt < this->N; ipt++)
    {
		if(!this->processed[ipt]) this->ExpandClusterOrder(ipt);
    }
	// Order chromosome and their reachDist in OPTICS
	// std::vector< std::pair< std::pair< chromosome*, int> , float > > OPTICS(this->N);
	// OPTICS.reserve(this->N);

    std::priority_queue< Pose, std::vector<Pose>, PoseClassifier::PoseClassifier > OPTICS;
    std::vector<Pose> vec;

	for(int i = 0; i < this->N; ++i)
	{
		// NEW
		// 
		if((this->points[i]).first != NULL && boost::math::isfinite(this->reachDist[i]) && this->order[i] <= this->N && this->order[i] >= 0)
		{
			Pose::Pose Pose((this->points[i]).first, i, this->order[i], this->reachDist[i], this->Population->Temperature);
			OPTICS.push(Pose);
            vec.push_back(Pose);
		}
		// OLD
		// 
		// OPTICS pairs contain :
		//   first  -> pair<chromosome*, index>
		//   second -> float reachDist
        // *** RE-DO THIS WHOLE FUCKING INSERTION SECTION ***
        // VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
        // place vector<>::iterator it @ the order emplacement in OPTICS
        // std::vector< std::pair< std::pair< chromosome*, int >, float> >::iterator it = OPTICS.begin()+this->order[i];
        
        // std::pair< chromosome*, int > newIntPair = std::make_pair((this->points[i]).first, i);
        // std::pair< std::pair< chromosome*, int >, float> newPair = std::make_pair( newIntPair, this->reachDist[i]);
        // if(newIntPair.first != NULL && (newIntPair.first)->genes != NULL && newIntPair.second >= 0 && newIntPair.second < this->N)
        // {
        //     it = OPTICS.erase(it);
        //     it = OPTICS.insert(it, newPair);
        // }
	}
    std::sort(vec.begin(),vec.end(),PoseClassifier::PoseClassifier());

	// Build BindingModes
    int i = 0; // used to have an idea of the number of loop completed (iterators are less convenient for that information while debugging)
	// NEW
//	while(!OPTICS.empty())
//	{
//		BindingMode::BindingMode current(this->Population);
//		while((OPTICS.top()).reachDist < 0.4 && !OPTICS.empty())
//		{
//			Pose currPose = OPTICS.top();
//			OPTICS.pop(); ++i;
//			current.add_Pose(currPose);
//		}
//		while((OPTICS.top()).reachDist >= 0.4 && !OPTICS.empty())
//		{
//			OPTICS.pop(); ++i;
//		}
//		if(current.get_BindingMode_size() >= 2)
//		{
//			this->Population->add_BindingMode(current);
//		}
//	}
	// OLD
 	BindingMode::BindingMode current(this->Population);
	for(std::vector< Pose >::iterator it = vec.begin(); it != vec.end(); ++i, ++it)
	{
        if(it->reachDist < 0.4) current.add_Pose(*it);
	   	if(it->reachDist >= 0.4)
        {
			if(current.get_BindingMode_size() >= 2)
			{
				this->Population->add_BindingMode(current);
			}
            current.clear_Poses();
        }
	}
}
std::vector<float> FastOPTICS::Vectorized_Chromosome(chromosome* chrom)
{
	std::vector<float> vChrom(this->nDimensions, 0.0f);
	// getting nDim-2 because the Dim=0 fills 3 memory cases
	for(int j = 0; j < this->nDimensions-2; ++j)
	{
		if(j == 0) //  building the first 3 comp. from genes[0] which are CartCoord x,y,z
		{
			for(int i = 0; i < 3; ++i)
			{
				vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].coor[i] - this->FA->ori[i]);
//                vChrom[i] = static_cast<float>( this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].coor[i] );
			}
		}
		else
		{
			// j+2 is used from {j = 1 to N} to build further comp. of genes[j]
			// vChrom[j+2] = static_cast<float>(genetoic(&gene_lim[j], (*chrom).genes[j].to_int32));
			vChrom[j+2] = static_cast<float>(RandomDouble( (*chrom).genes[j].to_int32 ));
		}
	}
	return vChrom;
}

void FastOPTICS::ExpandClusterOrder(int ipt)
{
    std::priority_queue< ClusterOrdering, std::vector<ClusterOrdering>, ClusterOrderingComparator::ClusterOrderingComparator > queue;
	ClusterOrdering tmp(ipt,0,1e6f);
	queue.push(tmp);

    while(!queue.empty() /*&& FastOPTICS::iOrder < this->N*/)
	{
		ClusterOrdering current = queue.top();
		queue.pop();
		int currPt = current.objectID;
		this->order[FastOPTICS::iOrder] = currPt;
		
		if(this->processed[currPt] == true) continue;
		// incrementing STATIC rank ordering ()
		FastOPTICS::iOrder++;
		this->processed[currPt] = true;
		
		float coredist = this->inverseDensities[currPt];
		for( std::vector<int>::iterator it = this->neighbors[currPt].begin(); it != this->neighbors[currPt].end(); ++it)
		{
			int iNeigh = *it;
			if(this->processed[iNeigh] == true) continue;

			float nrdist = this->compute_distance(points[iNeigh], points[currPt]);

			if(coredist > nrdist)
				nrdist = coredist;
			if(isUndefinedDist(this->reachDist[iNeigh]))
				this->reachDist[iNeigh] = nrdist;
			else if(nrdist < this->reachDist[iNeigh])
				this->reachDist[iNeigh] = nrdist;
			tmp = ClusterOrdering::ClusterOrdering(iNeigh, currPt, nrdist);
			queue.push(tmp);
		}
	}
}

float FastOPTICS::compute_distance(std::pair< chromosome*,std::vector<float> > & a, std::pair< chromosome*,std::vector<float> > & b)
{
	float distance = 0.0f;
	// chromosome* aChrom = a.first;
	// chromosome* branch = b.first;
	// std::vector<float> aVec(a.second);
	// std::vector<float> bVec(b.second);

	// insert distance calculation below
	for(int i = 0; i < this->nDimensions; ++i)
	{
		// distance += (aVec[i]-bVec[i])*(bVec[i]-aVec[i]);
		// x.second gives the vector reference of this->nDimensions size
        // float raw_distance = (a.second[i]-b.second[i]) * (b.second[i]-a.second[i]);
        // if(boost::math::isinf(raw_distance)) continue;
        // else if(boost::math::isnan(raw_distance)) continue;
        // else if(boost::math::isfinite(raw_distance)) distance += fabs(raw_distance);

        // simple
        // distance = (a.second[i]-b.second[i]) * (b.second[i]-a.second[i]);
        distance += (a.second[i]-b.second[i]) * (a.second[i]-b.second[i]);
	}
   	if(boost::math::isfinite(distance)) return sqrtf(distance);
	else return UNDEFINED_DIST;
}

/*****************************************\
			RandomProjections
\*****************************************/

// STATIC variables declaration
int const RandomProjectedNeighborsAndDensities::logOProjectionConstant;
float RandomProjectedNeighborsAndDensities::sizeTolerance;

// Constructor
RandomProjectedNeighborsAndDensities::RandomProjectedNeighborsAndDensities(std::vector< std::pair< chromosome*,std::vector<float> > >& inPoints, int minSplitSize, FastOPTICS* top)
{
	this->top = top;
	this->minSplitSize = minSplitSize;
	RandomProjectedNeighborsAndDensities::sizeTolerance = static_cast<float>(2.0f/3.0f);

	if( inPoints.size() < 1 )
	{
		this->N = 0;
		this->nDimensions = 0;
		this->nProject1D = 0;
		return;
	}
	else
	{
        this->points = inPoints;
		this->N = this->points.size();
		this->nDimensions = this->top->nDimensions;
		this->nPointsSetSplits = static_cast<int>(RandomProjectedNeighborsAndDensities::logOProjectionConstant * log(this->N * this->nDimensions + 1)/log(2));
		this->nProject1D = static_cast<int>(RandomProjectedNeighborsAndDensities::logOProjectionConstant * log(this->N * this->nDimensions + 1)/log(2));
		// line below calls copy-constructor
	}
	this->projectedPoints.reserve(this->nProject1D);
	for(int i = 0; i < this->nProject1D; ++i)
	{
        this->projectedPoints.push_back(std::vector<float>(this->N));
	}
//	this->splitsets.reserve(this->nPointsSetSplits);
//	for(int j = 0; j< this->nPointsSetSplits; ++j)
//	{
		// this->splitsets.push_back(std::vector<int>());
//	}
}

void RandomProjectedNeighborsAndDensities::computeSetBounds(std::vector< int > & ptList)
{
	std::vector< std::vector<float> > tempProj(this->nProject1D);
	// perform projection of points
	for(int j = 0; j<this->nProject1D; ++j)
	{
		std::vector<float> currentRp = this->Randomized_Normalized_Vector();

		int k = 0;
		std::vector<int>::iterator it = ptList.begin();
		while(it != ptList.end())
		{
			float sum = 0.0f;
			std::vector<float>::iterator vecPt = this->points[(*it)].second.begin();
			std::vector<float>::iterator currPro = (this->projectedPoints[j]).begin();
			for(int m = 0; m < this->nDimensions; ++m)
			{
				sum += currentRp[m] * vecPt[m];
			}
			currPro[k] = sum;
			++k;
            ++it;
		}
	}

	// Split Points Set
	std::vector<int> projInd(this->nProject1D);
	projInd.reserve(this->nProject1D);
	for(int j = 0; j < this->nProject1D; ++j) 
	{
//		projInd.push_back(j);
        projInd[j] = j;
	}
	for(int avgP = 0; avgP < this->nPointsSetSplits; avgP++)
	{
		// Shuffle projections
		for(int i = 0; i < this->nProject1D; ++i)
        {
//			tempProj.push_back(this->projectedPoints[i]);
            tempProj[i] = this->projectedPoints[i];
        }
		std::random_shuffle(projInd.begin(), projInd.end());
		std::vector<int>::iterator it = projInd.begin();
		int i = 0;
		while(it != projInd.end())
		{
			int cind = (*it);
			// look this line to be sure that the vector is pushed in this->projectedPoints
			this->projectedPoints[cind] = tempProj[i];
            i++;
            it++;
		}

		//split points set
		int nPoints = ptList.size();
		std::vector<int> ind(nPoints);
		ind.reserve(nPoints);
		for(int l = 0; l < nPoints; ++l)
		{
			ind[l] = l;
		}
		this->SplitUpNoSort(ind,0);
	}
}

void RandomProjectedNeighborsAndDensities::SplitUpNoSort(std::vector< int >& ind, int dim)
{
	int nElements = ind.size();
	dim = dim % this->nProject1D;
	std::vector<float>::iterator tProj = (this->projectedPoints[dim]).begin();
	int splitPos;

	// save set such that used for density or neighborhood computation
	if(nElements > this->minSplitSize*(1 - RandomProjectedNeighborsAndDensities::sizeTolerance) && nElements < this->minSplitSize*(1 + RandomProjectedNeighborsAndDensities::sizeTolerance))
	{
		std::vector<float> cpro(nElements);
		for(int i = 0; i < nElements; ++i)	
		{
//			cpro.push_back(tProj[ind[i]]);
            cpro[i] = tProj[ind[i]];
		}
		// sprting cpro[] && ind[] concurrently
		this->quicksort_concurrent_Vectors(cpro, ind, 0, nElements-1);
		this->splitsets.push_back(ind);
	}

	// compute splitting element
	if(nElements > this->minSplitSize)
	{
		//pick random splitting element based on position
		int randInt = static_cast<int>(std::floor( (RandomDouble()*static_cast<double>(nElements)) ));
		float rs = tProj[ind[randInt]];
		int minInd = 0;
		int maxInd = nElements - 1;
		while(minInd < maxInd)
		{
			float currEle = tProj[ind[minInd]];
			if(currEle > rs)
			{
				while(minInd < maxInd && tProj[ind[maxInd]] > rs) maxInd--;
				if(minInd == maxInd) break;
				
				int currInd = ind[minInd];
				ind[minInd]=ind[maxInd];
				ind[maxInd]=currInd;
				maxInd--;
			}
			minInd++;
		}
		if( minInd == nElements-1 ) minInd = nElements/2;

		// split set recursively
		splitPos = minInd + 1;
		
		// std::vector<int> ind2(splitPos);
		std::vector<int> ind2(splitPos);
		for(int l = 0; l < splitPos; ++l)
		{
			 ind2[l] = ind[l];
//			ind2.push_back(ind[l]);
		}
		this->SplitUpNoSort(ind2,dim+1);
		
		// std::vector<int> ind3(nElements - splitPos);
        ind2 = std::vector<int>( nElements-splitPos );
		for(int l = 0; l < nElements-splitPos; ++l)
		{
			 ind2[l] = ind[l+splitPos];
//			ind3.push_back(ind[l+splitPos]);
		}
		this->SplitUpNoSort(ind2,dim+1);
	}
}

void RandomProjectedNeighborsAndDensities::getInverseDensities(std::vector< float > & inverseDensities)
{
	inverseDensities.reserve(this->N);
	std::vector<int> nDists(this->N);
	nDists.reserve(this->N);
//	for(std::vector< std::vector< int> >::iterator it1 = this->splitsets.begin(); it1 != this->splitsets.end(); ++it1)
    for(int i = 0; i < this->splitsets.size(); ++i)
	{
//		std::vector<int>::iterator pinSet = it1->begin();
        std::vector<int> & pinSet = this->splitsets.at(i);
		
//		int len = it1->size();
        int len = pinSet.size();
		int indoff = static_cast<int>(round(len/2));
		int oldind = pinSet[indoff];
		for(int i = 0; i < len; ++i)
		{
			int ind = pinSet[i];
			if(ind == indoff) continue;
			// CHOOSE BETWEEN OF THE 2 CALLS BELOW FOR compute_distance() IMPLEMENTATION
			float dist = this->top->compute_distance(this->points[ind],this->points[oldind]);
			inverseDensities[oldind] += dist;
			nDists[oldind]++;
			inverseDensities[ind] += dist;
			nDists[ind]++;
		}
	}
	for(int l = 0; l < this->N; ++l)
	{
		if(nDists[l] == 0) inverseDensities[l] = UNDEFINED_DIST;
		else inverseDensities[l] /= nDists[l];
	}
}

void RandomProjectedNeighborsAndDensities::getNeighbors(std::vector< std::vector< int > > & neighs)
{
	neighs.reserve(this->N);
	for(int l = 0; l < this->N; l++)
	{
		std::vector<int> list;
		// std::vector<int> list(this->N,0);
		// list.reserve(this->N);
		neighs.push_back(list);
	}

	// go through all sets
	for(std::vector< std::vector< int > >::iterator it1 = this->splitsets.begin(); it1 != this->splitsets.end(); ++it1)
	{
		// for each set (each projected line)
		std::vector<int>::iterator pinSet = it1->begin();
		int len = it1->size();
		int ind = pinSet[0];
		int indoff = static_cast<int>(round(len/2));
		int oldind = pinSet[indoff];

		// add all points as neighbors to middle point
		//  +
		// add the middle point to all other points in set
		for(int i = 0; i < len; ++i)
		{
			ind = pinSet[i];

			// The following block of code uses an iterator to check 
			std::vector<int> & cneighs = neighs.at(ind);
			// std::vector<int>::iterator itPos = std::find(cneighs.begin(),cneighs.end(),oldind);
			// if(itPos == cneighs.end()) //element not found in cneigh
			if( !std::binary_search(cneighs.begin(), cneighs.end(), oldind) )
			{
				// if(*itPos > cneighs.size())
				// 	cneighs.push_back(oldind);
				// else
				// 	cneighs.insert(itPos, oldind);
				cneighs.push_back(oldind);
				std::sort(cneighs.begin(),cneighs.end());
				cneighs.erase(std::unique(cneighs.begin(), cneighs.end()), cneighs.end());
			}

			std::vector<int> & cneighs2 = neighs.at(oldind);
			// itPos = std::find(cneighs2.begin(), cneighs2.end(), ind);
			// if(itPos == cneighs2.end()) // element not found in cneighs2
			if( !std::binary_search(cneighs2.begin(), cneighs2.end(), ind) )
			{
				// if(*itPos > cneighs2.size())
				// 	cneighs2.push_back(ind);
				// else
				// 	cneighs2.insert(itPos, ind);
				cneighs2.push_back(oldind);
				std::sort(cneighs2.begin(),cneighs2.end());
				cneighs2.erase(std::unique(cneighs2.begin(), cneighs2.end()), cneighs2.end());
			}
		}

	}
}

std::vector<float> RandomProjectedNeighborsAndDensities::Randomized_Normalized_Vector()
{
	std::vector<float> vChrom(this->nDimensions, 0.0f);
	
	float sum = 0.0f;
	for(int j = 0; j < this->nDimensions-2; ++j)
	{
		int intGene = this->Dice();
		
		double doubleGene = RandomDouble(intGene);
		if(j == 0) //  building the first 3 comp. from genes[0] which are CartCoord x,y,z
		{
            double doubleGeneIC = genetoic(&this->top->gene_lim[j],intGene);
			for(int i = 0; i < 3; ++i)
			{
				vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].coor[i] - this->top->FA->ori[i]);
				sum += vChrom[i]*vChrom[i];
			}
		}
		else
		{
			// j+2 is used from {j = 1 to N} to build further comp. of genes[j]
			vChrom[j+2] = static_cast<float>(doubleGene);
			sum += vChrom[j+2]*vChrom[j+2];
		}
	}
	sum = sqrtf(sum);
	for(int k = 0; k < this->nDimensions; ++k) { vChrom[k]/=sum; }

	return vChrom;
}

void RandomProjectedNeighborsAndDensities::quicksort_concurrent_Vectors(std::vector<float>& data, std::vector<int>& index, int beg, int end)
{
	int l, r, p;
	float pivot;
	std::vector<float>::iterator xData, yData;
	std::vector<int>::iterator xIndex, yIndex;
	while(beg < end)
	{
		l = beg; p = beg + (end-beg)/2; r = end;
		pivot = data[p];
		
		while(1)
		{
			while( (l<=r) && QS_ASC(data[l],pivot) <= 0 ) ++l;
			while( (l<=r) && QS_ASC(data[r],pivot)  > 0 ) --r;
	
			if (l > r) break;
			xData = data.begin()+l; yData = data.begin()+r;
			xIndex = index.begin()+l; yIndex = index.begin()+r;
			this->swap_element_in_vectors(xData, yData, xIndex, yIndex);
			if (p == r) p = l;
			++l; --r;
		}
		xData = data.begin()+p; yData = data.begin()+r;
		xIndex = index.begin()+p; yIndex = index.begin()+r;
		this->swap_element_in_vectors(xData, yData, xIndex, yIndex);

		--r;

		if( (r-beg) < (end-l) )
		{
			this->quicksort_concurrent_Vectors(data, index, beg, r);
			beg = l;
		}
		else
		{
			this->quicksort_concurrent_Vectors(data, index, l, end);
			end = r;
		}
	}
}

void RandomProjectedNeighborsAndDensities::swap_element_in_vectors(std::vector<float>::iterator xData, std::vector<float>::iterator yData, std::vector<int>::iterator xIndex, std::vector<int>::iterator yIndex)
{
	float tData = *xData; *xData = *yData; *yData = tData;
	int tIndex = *xIndex; *xIndex = *yIndex; *yIndex = tIndex;
}

// This function generates a RandomInt32 who can be used as *genes->to_int32 value
int RandomProjectedNeighborsAndDensities::Dice()
{
	unsigned int tt = static_cast<unsigned int>(time(0));
	srand(tt);
	RNGType rng(tt);
	boost::uniform_int<> one_to_max_int32( 0, MAX_RANDOM_VALUE );
	boost::variate_generator< RNGType, boost::uniform_int<> > dice(rng, one_to_max_int32);
	return dice();
}