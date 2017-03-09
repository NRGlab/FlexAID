#include "ColonyEnergy.h"
#include "gaboom.h"

/*****************************************\
				ColonyEnergy
\*****************************************/
// Constructor and Algorithm main+only call
ColonyEnergy::ColonyEnergy(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, BindingPopulation& Population, int nPoints)
{	
// Declarations
	///////////////////////////////////////////////////////
	// Entropy
	this->Population = &Population;
	// FlexAID
    this->FA = FA;
    this->GB = GB;
    this->VC = VC;
    this->cleftgrid = cleftgrid;
    this->atoms = atoms;
    this->residue = residue;
    this->chroms = chrom;
    this->gene_lim = gen_lim;
	this->N = static_cast<int>( this->Population->Poses.size() );
	this->dist_threshold = this->FA->cluster_rmsd;
	// this->dist_threshold = this->FA->cluster_rmsd*(2 - RandomProjectedNeighborsColonyEnergy::sizeTolerance);
	// ColonyEnergy
    this->nDimensions = this->FA->num_het_atm*3;	// use with RandomProjectedNeighborsColonyEnergy()
    // this->nDimensions = this->FA->npar + 2; 	// use with Vectorized_Chromosome()
    
    this->minPoints = nPoints;
    
	this->processed.reserve(this->N);
	this->colonyEnergy.reserve(this->N);
	this->points.reserve(this->N);
	this->neighbors.reserve(this->N);
    
    for(std::vector<Pose>::iterator iPose = this->Population->Poses.begin(); iPose != this->Population->Poses.end(); ++iPose)
	{
        // default values initialization
        this->processed.push_back(false);
        double Pi = ( iPose->boltzmann_weight / this->Population->PartitionFunction );
        double CFdS = Pi*iPose->CF - ( this->Population->Temperature * (-1 * Pi * log(Pi)) );
        this->colonyEnergy.push_back( CFdS );
        // std::pair<first, second> is pushed to this->points[]
        // 	first  -> chromosome* pChrom (pointer to chromosome)
        // 	second -> vChrom[FA->npar+n] (vectorized chromosome)
        this->points.push_back( std::make_pair( iPose->chrom, iPose->vPose) ) ;
	}
}

void ColonyEnergy::Execute_ColonyEnergy(char* end_strfile, char* tmp_end_strfile)
{
	std::vector< int > ptInd;
    ptInd.reserve(this->Population->Poses.size());

    for(int k = 0; k < this->Population->Poses.size(); ++k)
    {
        ptInd.push_back(k);
    }

	RandomProjectedNeighborsColonyEnergy MultiPartition = RandomProjectedNeighborsColonyEnergy(this->points, this->minPoints, this);
	
	// Fast— section of FastOPTICS to find neighbors
	MultiPartition.computeSetBounds(ptInd);
    
    // getNeighbors also updates CFdS attribute of the Poses in Population
	MultiPartition.getNeighbors(this->neighbors);

	// sort Poses by CFdS
    std::sort(this->Population->Poses.begin(), this->Population->Poses.end(), PoseRanker());
}


/*****************************************\
			RandomProjections
\*****************************************/

// STATIC variables declaration
int const RandomProjectedNeighborsColonyEnergy::logOProjectionConstant = 20;
float RandomProjectedNeighborsColonyEnergy::sizeTolerance = static_cast<float>(2.0f/3.0f);

// Constructor
RandomProjectedNeighborsColonyEnergy::RandomProjectedNeighborsColonyEnergy(std::vector< std::pair< chromosome*,std::vector<float> > >& inPoints, int minSplitSize, ColonyEnergy* top)
{
	this->top = top;
	this->minSplitSize = minSplitSize;

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
		this->N = static_cast<int>(this->points.size());
		this->nDimensions = this->top->nDimensions;
		
		// this->nPointsSetSplits = static_cast<int>(RandomProjectedNeighborsColonyEnergy::logOProjectionConstant * log(static_cast<float>(this->N * (this->top->FA->npar+2) + 1))/log(2.f));
		// this->nProject1D = static_cast<int>(RandomProjectedNeighborsColonyEnergy::logOProjectionConstant * log(static_cast<float>(this->N * (this->top->FA->npar+2) + 1))/log(2.f));
		
		this->nPointsSetSplits = static_cast<int>(RandomProjectedNeighborsColonyEnergy::logOProjectionConstant * log(static_cast<float>(this->N * (this->nDimensions) + 1))/log(2.f));
		this->nProject1D = static_cast<int>(RandomProjectedNeighborsColonyEnergy::logOProjectionConstant * log(static_cast<float>(this->N * (this->nDimensions) + 1))/log(2.f));
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

void RandomProjectedNeighborsColonyEnergy::computeSetBounds(std::vector< int > & ptList)
{
	std::vector< std::vector<float> > tempProj(this->nProject1D);
	// perform projection of points
	for(int j = 0; j < this->nProject1D; ++j)
	{
        // std::vector<float> currentRp( this->Randomized_InternalCoord_Vector() );
        std::vector<float> currentRp( this->Randomized_CartesianCoord_Vector() );
        // std::vector<float> currentRp( this->Randomly_Selected_Chromosome() );
        
		int k = 0;
		std::vector<int>::iterator it = ptList.begin();
		while(it != ptList.end())
		{
			float sum = 0.0f;
			// std::vector<float>::iterator vecPt = this->points[(*it)].second.begin();
			std::vector<float> vecPt(this->points[(*it)].second);
			std::vector<float>::iterator currPro = (this->projectedPoints[j]).begin();
			
			for(int m = 0; m < this->nDimensions; ++m)
				sum += currentRp[m] * vecPt[m];

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
//        std::random_shuffle(projInd.begin(), projInd.end(), ([](int n) { return rand() % n; }) );
        std::random_shuffle(projInd.begin(), projInd.end());
		
		int i = 0;
        for(std::vector<int>::iterator it = projInd.begin(); it != projInd.end(); ++it, i++)
		{
			int cind = (*it);
			// look this line to be sure that the vector is pushed in this->projectedPoints
			this->projectedPoints[cind] = tempProj[i];
		}

		//split points set
		int nPoints = static_cast<int>(ptList.size());
		std::vector<int> ind(nPoints);
		ind.reserve(nPoints);
		for(int l = 0; l < nPoints; ++l)
		{
			ind[l] = l;
		}
		this->SplitUpNoSort(ind,0);
	}
}

void RandomProjectedNeighborsColonyEnergy::SplitUpNoSort(std::vector< int >& ind, int dim)
{
	int nElements = static_cast<int>(ind.size());
	// dim = dim % this->nProject1D;
	dim = rand() % this->nProject1D;
	std::vector<float>::iterator tProj = (this->projectedPoints[dim]).begin();
	int splitPos;

	// save set such that used for density or neighborhood computation
	if(nElements > this->minSplitSize*(1 - RandomProjectedNeighborsColonyEnergy::sizeTolerance) && nElements < this->minSplitSize*(1 + RandomProjectedNeighborsColonyEnergy::sizeTolerance))
	{
		std::vector<float> cpro(nElements);
		for(int i = 0; i < nElements; ++i)	
		{
//			cpro.push_back(tProj[ind[i]]);
            cpro[i] = tProj[ind[i]];
		}
		// sprting cpro[] && ind[] concurrently
		this->quicksort_concurrent_Vectors(cpro, ind, 0, nElements-1);
		
		// push back the split points set
		this->splitsets.push_back(ind);
	}

	// compute splitting element
	if(nElements > this->minSplitSize)
	{
		//pick random splitting element based on position
		// int randInt = static_cast<int>(std::floor( (RandomDouble()*static_cast<double>(nElements)) ));
		int randInt = rand() % nElements;
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
				ind[minInd] = ind[maxInd];
				ind[maxInd] = currInd;
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
		
        std::vector<int> ind3(nElements - splitPos);
//        ind2 = std::vector<int>( nElements-splitPos );
		for(int l = 0; l < (nElements-splitPos); ++l)
		{
			 ind3[l] = ind[l+splitPos];
//			ind3.push_back(ind[l+splitPos]);
		}
		this->SplitUpNoSort(ind3,dim+1);
	}
}

void RandomProjectedNeighborsColonyEnergy::getNeighbors(std::vector< std::vector< int > > & neighs)
{
	neighs.reserve(this->N);
	for(int l = 0; l < this->N; l++)
	{
		std::vector<int> list;
		neighs.push_back(list);
	}

	// go through all sets
	for(std::vector< std::vector< int > >::iterator it1 = this->splitsets.begin(); it1 != this->splitsets.end(); ++it1)
	{
		// for each set (each projected line)
		std::vector<int>::iterator pinSet = it1->begin();
		int len = static_cast<int>(it1->size());
		int ind = pinSet[0];
		int indoff = static_cast<int>(round_it(len/2));
		int oldind = pinSet[indoff];

		// add all points as neighbors to middle point
		//  +
		// add the middle point to all other points in set
		for(int i = 0; i < len; ++i)
		{
			ind = pinSet[i];

			// check if the distance is <= to dist_threshold before condiering a point as a neighbor
			if( this->top->compute_distance(this->top->points[ind], this->top->points[oldind]) > this->top->dist_threshold ) continue;
			
			// The following block of code uses an iterator to add all points as neighbors to the middle point
			std::vector<int> & cneighs = neighs.at(ind);
			
			if( !std::binary_search(cneighs.begin(), cneighs.end(), oldind) ) //only add point if not a neighbor already
			{
                this->top->Population->Poses[ind].CFdS += this->top->colonyEnergy[oldind];
                
				cneighs.push_back(oldind);
				std::sort(cneighs.begin(),cneighs.end());
//              cneighs.erase(std::unique(cneighs.begin(), cneighs.end()), cneighs.end());
                
			}
			// The following block of code adds the middle point as neighbor to all other points in set
			std::vector<int> & cneighs2 = neighs.at(oldind);
			if( !std::binary_search(cneighs2.begin(), cneighs2.end(), ind) ) //only add point if not a neighbor already
			{
                this->top->Population->Poses[oldind].CFdS += this->top->colonyEnergy[ind];
                
				cneighs2.push_back(ind);
				std::sort(cneighs2.begin(),cneighs2.end());
//				cneighs2.erase(std::unique(cneighs2.begin(), cneighs2.end()), cneighs2.end());
			}
		}

	}
}

std::vector<float> ColonyEnergy::Vectorized_Chromosome(chromosome* chrom)
{
    float norm = 0.0f;
	std::vector<float> vChrom(this->nDimensions, 0.0f);
	// getting nDim-2 because the Dim=0 fills 3 memory cases
	for(int j = 0; j < this->nDimensions-2; ++j)
	{
		if(j == 0) //  building the first 3 comp. from genes[0] which are CartCoord x,y,z
		{
			for(int i = 0; i < 3; ++i)
			{
				// vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].coor[i] - this->FA->ori[i]);
				if(i == 0)
				{
                    // vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].coor[i] - this->FA->ori[i]);
					vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].dis);
					// vChrom[i] *= vChrom[i];
				}
				if(i == 1)
				{
					// vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].coor[i] - this->FA->ori[i]);
					vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].ang);
					// vChrom[i] = static_cast<float>( RandomDouble( (*chrom).genes[j].to_int32) );
					// vChrom[i] *= vChrom[i];
				}
				if(i == 2)
				{
                    // vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].coor[i] - this->FA->ori[i]);
					vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].dih);
					// vChrom[i] = static_cast<float>( genetoic(&this->gene_lim[i],(*chrom).genes[j].to_int32) );
					// vChrom[i] *= vChrom[i];
				}
                norm += vChrom[i]*vChrom[i];
			}
		}
		else
		{
			// j+2 is used from {j = 1 to N} to build further comp. of genes[j]
			// vChrom[j+2] = static_cast<float>(genetoic(&gene_lim[j], (*chrom).genes[j].to_int32));
			vChrom[j+2] = static_cast<float>((*chrom).genes[j].to_ic);
			// vChrom[j+2] = static_cast<float>( RandomDouble( (*chrom).genes[j].to_int32) );
            norm += vChrom[j+2]*vChrom[j+2];
		}
	}
    
  // norm = sqrtf(norm);
  // for(int k = 0; k < this->nDimensions; ++k) { vChrom[k]/=norm; }
   
   return vChrom;
}

std::vector<float> ColonyEnergy::Vectorized_Cartesian_Coordinates(int chrom_index)
{
	int i = 0,j = 0,l = 0,m = 0;
	int cat;
	int rot;

	uint grd_idx;
	int normalmode=-1;
	int rot_idx=0;

    std::vector<float> vChrom(this->nDimensions);

	int npar = this->GB->num_genes;
	
	j = chrom_index;

	for(i=0;i<npar;i++){ this->FA->opt_par[i] = this->chroms[j].genes[i].to_ic; }

	for(i=0;i<npar;i++)
	{
		//printf("[%8.3f]",FA->opt_par[i]);
  
		if(this->FA->map_par[i].typ==-1) 
		{ //by index
			grd_idx = (uint)this->FA->opt_par[i];
			//printf("this->FA->opt_par(index): %d\n", grd_idx);
			//PAUSE;
			this->atoms[this->FA->map_par[i].atm].dis = this->cleftgrid[grd_idx].dis;
			this->atoms[this->FA->map_par[i].atm].ang = this->cleftgrid[grd_idx].ang;
			this->atoms[this->FA->map_par[i].atm].dih = this->cleftgrid[grd_idx].dih;

		}
		else if(this->FA->map_par[i].typ == 0)
		{
			this->atoms[this->FA->map_par[i].atm].dis = (float)this->FA->opt_par[i];
		}
		else if(this->FA->map_par[i].typ == 1)
		{
			this->atoms[this->FA->map_par[i].atm].ang = (float)this->FA->opt_par[i];
		}
		else if(this->FA->map_par[i].typ == 2)
		{
			this->atoms[this->FA->map_par[i].atm].dih = (float)this->FA->opt_par[i];

			j=this->FA->map_par[i].atm;
			cat=this->atoms[j].rec[3];
			if(cat != 0)
			{
				while(cat != this->FA->map_par[i].atm)
				{
					this->atoms[cat].dih=this->atoms[j].dih + this->atoms[cat].shift; 
					j=cat;
					cat=this->atoms[j].rec[3];
				}
			}
		}else if(this->FA->map_par[i].typ == 3)
		{ //by index
			grd_idx = (uint)this->FA->opt_par[i];

			// serves as flag , but also as grid index
			normalmode=grd_idx;

		}else if(this->FA->map_par[i].typ == 4)
		{
			rot_idx = (int)(this->FA->opt_par[i]+0.5);

			this->residue[this->atoms[this->FA->map_par[i].atm].ofres].rot=rot_idx;
		}
  
	}

	if(normalmode > -1) alter_mode(this->atoms,this->residue,this->FA->normal_grid[normalmode],this->FA->res_cnt,this->FA->normal_modes);

	/* rebuild cartesian coordinates of optimized residues*/
    for(i=0;i<this->FA->nors;i++) buildcc(this->FA,this->atoms,this->FA->nmov[i],this->FA->mov[i]);

	// residue that is optimized geometrically (ligand)
	l=this->atoms[this->FA->map_par[0].atm].ofres;

	rot=this->residue[l].rot;
    m=0;
	for(i=this->residue[l].fatm[rot];i<=this->residue[l].latm[rot];i++)
	{
		for(j=0;j<3;j++) vChrom[m*3+j] = this->atoms[i].coor[j];
        ++m;
	}
	return vChrom;
}

std::vector<float> RandomProjectedNeighborsColonyEnergy::Randomized_InternalCoord_Vector()
{
    std::vector<float> vChrom(this->nDimensions);

	float sum = 0.0f;
	for(int j = 0; j < this->nDimensions-2; ++j)
	{
		if(j == 0) //  building the first 3 comp. from genes[0] which are CartCoord x,y,z
		{
            double doubleGeneIC = genetoic(&this->top->gene_lim[j],roll_die());
            for(int i = 0; i < 3; ++i)
			{
				// vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].coor[i] - this->top-÷>FA->ori[i]);
				if(i == 0)
				{
					// vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].coor[i] - this->top->FA->ori[i]);
					vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].dis);
//                    vChrom[i] *= vChrom[i];
				}
				if(i == 1)
				{
                    // vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].coor[i] - this->top->FA->ori[i]);
					vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].ang);
//                    vChrom[i] *= vChrom[i];
				}
				if(i == 2)
				{
                    // vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].coor[i] - this->top->FA->ori[i]);
					vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].dih);
//                    vChrom[i] *= vChrom[i];
				}
				sum += vChrom[i]*vChrom[i];
			}
		}
		else
		{
			// j+2 is used from {j = 1 to N} to build further comp. of genes[j]
			// vChrom[j+2] = static_cast<float>(RandomDouble(random_dice()));
			vChrom[j+2] = static_cast<float>(genetoic(&this->top->gene_lim[j],roll_die()));
			sum += vChrom[j+2]*vChrom[j+2];
		}
	}
	sum = sqrtf(sum);
	// for(int k = 0; k < this->nDimensions; ++k) { vChrom[k]/=sum; }

    
	return vChrom;
}

std::vector<float> RandomProjectedNeighborsColonyEnergy::Randomized_CartesianCoord_Vector()
{
   // float norm = 0.0f;

    int i = 0,j = 0,l = 0,m = 0;
	int cat;
	int rot;

	uint grd_idx;
	int normalmode=-1;
	int rot_idx=0;

    std::vector<float> vChrom(this->nDimensions);
    
    // set the numper of parameters
	int npar = this->top->GB->num_genes;

    // generate random ic for each DoF
	for(i=0;i<npar;i++){ this->top->FA->opt_par[i] = genetoic(&this->top->gene_lim[i],roll_die()); }

    // copy randomly generated parameters above to FA->atoms structure prior to builc cartesian coordinates (cc) for the randomized individual
	for(i=0;i<npar;i++)
	{
  
		if(this->top->FA->map_par[i].typ==-1) 
		{ //by index
			grd_idx = (uint)this->top->FA->opt_par[i];

			this->top->atoms[this->top->FA->map_par[i].atm].dis = this->top->cleftgrid[grd_idx].dis;
			this->top->atoms[this->top->FA->map_par[i].atm].ang = this->top->cleftgrid[grd_idx].ang;
			this->top->atoms[this->top->FA->map_par[i].atm].dih = this->top->cleftgrid[grd_idx].dih;

		}
		else if(this->top->FA->map_par[i].typ == 0)
		{
			this->top->atoms[this->top->FA->map_par[i].atm].dis = (float)this->top->FA->opt_par[i];
		}
		else if(this->top->FA->map_par[i].typ == 1)
		{
			this->top->atoms[this->top->FA->map_par[i].atm].ang = (float)this->top->FA->opt_par[i];
		}
		else if(this->top->FA->map_par[i].typ == 2)
		{
			this->top->atoms[this->top->FA->map_par[i].atm].dih = (float)this->top->FA->opt_par[i];

			j=this->top->FA->map_par[i].atm;
			cat=this->top->atoms[j].rec[3];
			if(cat != 0)
			{
				while(cat != this->top->FA->map_par[i].atm)
				{
					this->top->atoms[cat].dih=this->top->atoms[j].dih + this->top->atoms[cat].shift; 
					j=cat;
					cat=this->top->atoms[j].rec[3];
				}
			}
		}else if(this->top->FA->map_par[i].typ == 3)
		{ //by index
			grd_idx = (uint)this->top->FA->opt_par[i];

			// serves as flag , but also as grid index
			normalmode=grd_idx;

		}else if(this->top->FA->map_par[i].typ == 4)
		{
			rot_idx = (int)(this->top->FA->opt_par[i]+0.5);

			this->top->residue[this->top->atoms[this->top->FA->map_par[i].atm].ofres].rot=rot_idx;
		}
  
	}
    
    // process normal modes (if used in the simulation)
	if(normalmode > -1) alter_mode(this->top->atoms, this->top->residue, this->top->FA->normal_grid[normalmode], this->top->FA->res_cnt, this->top->FA->normal_modes);

	/* rebuild cartesian coordinates of optimized residues*/
    for(i=0;i<this->top->FA->nors;i++) buildcc(this->top->FA,this->top->atoms,this->top->FA->nmov[i],this->top->FA->mov[i]);

	// residue that is optimized geometrically (ligand)
	l=this->top->atoms[this->top->FA->map_par[0].atm].ofres;

	rot=this->top->residue[l].rot;
    m=0;
	for(i = this->top->residue[l].fatm[rot]; i <= this->top->residue[l].latm[rot]; i++)
	{
		for(j=0;j<3;j++)
		{
			vChrom[m*3+j] = this->top->atoms[i].coor[j];
            // norm += vChrom[m*3+j]*vChrom[m*3+j];
		}
        ++m;
	}
   //  norm = sqrtf(norm);
  	// for(i = 0; i < this->nDimensions; ++i) vChrom[i] /= norm;
	return vChrom;
}

int ColonyEnergy::get_minPoints() { return this->minPoints; }

std::vector<int> ColonyEnergy::get_neighbors_for_chrom(int chrom_index)
{
	std::vector<int> neighs(this->neighbors[chrom_index]);
	return neighs;
}
/*****************************************\

			RandomProjections
\*****************************************/

float ColonyEnergy::compute_distance(std::pair< chromosome*,std::vector<float> > & a, std::pair< chromosome*,std::vector<float> > & b)
{
	float distance = 0.0f;

	// simple distance calculation below
	for(int i = 0; i < this->nDimensions; ++i)
	{
		float tempDist = (a.second[i]-b.second[i]);
        distance +=  tempDist * tempDist;
	}

	return sqrtf(distance / static_cast<float>(this->FA->num_het_atm));
}

float ColonyEnergy::compute_vect_distance(std::vector<float> a, std::vector<float> b)
{
	float distance = 0.0f;

	// simple distance calculation below
	for(int i = 0; i < this->nDimensions; ++i)
    {
		float tempDist = (a[i]-b[i]);
        distance +=  tempDist * tempDist;
    }

	return sqrtf(distance / static_cast<float>(this->FA->num_het_atm));
}

std::vector<float> RandomProjectedNeighborsColonyEnergy::Randomly_Selected_Chromosome()
{
	int ind = rand() % this->N;
    std::vector<float> vChrom(this->top->Vectorized_Cartesian_Coordinates(ind));
	return vChrom;
}

void RandomProjectedNeighborsColonyEnergy::quicksort_concurrent_Vectors(std::vector<float>& data, std::vector<int>& index, int beg, int end)
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

void RandomProjectedNeighborsColonyEnergy::swap_element_in_vectors(std::vector<float>::iterator xData, std::vector<float>::iterator yData, std::vector<int>::iterator xIndex, std::vector<int>::iterator yIndex)
{
	float tData = *xData; *xData = *yData; *yData = tData;
	int tIndex = *xIndex; *xIndex = *yIndex; *yIndex = tIndex;
}

