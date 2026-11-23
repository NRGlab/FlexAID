#include "FOPTICS.h"
#include "gaboom.h"

#include <random>
#include <algorithm>
#include <chrono>

std::random_device rd;
std::mt19937 gen(rd());

struct RNG
{
    int operator() (int n) {
        return std::rand() / (1.0 + RAND_MAX) * n;
    }
};

bool definitelyGreaterThan(float a, float b, float epsilon)
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyLessThan(float a, float b, float epsilon)
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

//Normalizes any number to an arbitrary range
//by assuming the range wraps around when going below min or above max
float normalize_IC_interval(const genlim* gene_lim, float value)
{
	float start = gene_lim->min;
	float end = gene_lim->max;
	float width = end - start;
	float offsetValue = value - start; // value relative to 0
	return ( offsetValue - ( std::floor( offsetValue / width) * width) ) + start;
}

ClusterOrdering::ClusterOrdering(int id, int predID, float reach) : objectID(id), predecessorID(predID), reachability(reach)
{}

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
	if(this->reachability > rhs.reachability || isUndefinedDist(this->reachability))
			return false;
		else if(this->reachability < rhs.reachability)
			return true;
		if(this->objectID > rhs.objectID)
			return true;
		else if(this->objectID < rhs.objectID)
			return false;
		// if nothing else is true, return 0
		return 0;
}

inline bool const ClusterOrdering::operator> (const ClusterOrdering& rhs)
{
	if(this->reachability > rhs.reachability || isUndefinedDist(this->reachability))
		return true;
	else if(this->reachability < rhs.reachability)
			return false;
	if(this->objectID > rhs.objectID)
		return false;
	else if(this->objectID < rhs.objectID)
		return true;
		// if nothing else is true, return 0
		return 0;
}
/*****************************************\
				FastOPTICS
\*****************************************/
// Constructor and Algorithm main+only call
// int FastOPTICS::iOrder;
FastOPTICS::FastOPTICS(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, BindingPopulation& Population, int nPoints)
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
	this->N = nChrom;

	// FastOPTICS
    this->nDimensions = this->FA->num_het_atm*3;	// use with Vectorized_Cartesian_Coordinates()
    // this->nDimensions = this->FA->npar + 2; 	// use with Vectorized_Chromosome()

    this->minPoints = nPoints;

    // FastOPTICS::iOrder = 0;
    this->iOrder = 0;
	this->order.reserve(this->N);
	this->reachDist.reserve(this->N);
	this->processed.reserve(this->N);
	this->inverseDensities.reserve(this->N);
	this->points.reserve(this->N);
	this->neighbors.reserve(this->N);

	for(int i = 0; i < this->N; ++i)
	{
		// need to transform the chromosomes into vector f
		// std::vector<float> vChrom(this->FastOPTICS::Vectorized_Chromosome(&chrom[i])); // Copy constructor
		std::vector<float> vChrom( this->Vectorized_Cartesian_Coordinates(i) ); // Copy constructor
		if(vChrom.size() == this->nDimensions)
		{
			// std::pair<first, second> is pushed to this->points[]
			// 	first  -> chromosome* pChrom (pointer to chromosome)
			// 	second -> vChrom[FA->npar+n] (vectorized chromosome)
			this->order.push_back(-1);
			this->reachDist.push_back(UNDEFINED_DIST);
			this->processed.push_back(false);
			this->inverseDensities.push_back(0.0f);
			this->points.push_back( std::make_pair( (chromosome*)&chrom[i], vChrom) ) ;
		}
	}

};

void FastOPTICS::Execute_FastOPTICS(char* end_strfile, char* tmp_end_strfile)
{
	// vector of point indexes
	std::vector< int > ptInd;
    ptInd.reserve(this->N);
    for(int k = 0; k < this->N; ++k)
    {
        ptInd.push_back(k);
    }

	// Build object, compute projections, density estimates and density neighborhoods (in serial order of function calls below)
	RandomProjectedNeighborsAndDensities::RandomProjectedNeighborsAndDensities MultiPartition(this->points, this->minPoints, this); // use minPoints amount of random projections
	MultiPartition.computeSetBounds(ptInd);
	MultiPartition.getInverseDensities(this->inverseDensities);
	MultiPartition.getNeighbors(this->neighbors);

	// Compute OPTICS ordering
	for(int ipt = 0; ipt < this->N; ipt++) 		// starting from 0
//	for(int ipt = this->N-1; ipt >= 0; --ipt)	// starting from this->N-1
    {
		if(!this->processed[ipt]) this->ExpandClusterOrder(ipt);
    }

    // Would it be useful to normalize 'reachability distance' ?
    // this->normalizeDistances();
    // output the projected distances

    MultiPartition.output_projected_distance(end_strfile, tmp_end_strfile);

	// Order chromosome and their reachDist in OPTICS
    //  points pairs contain :
    //   first  -> pair<chromosome*, index>
    //   second -> float reachDist
	for(int i = 0; i < this->N; ++i)
	{
		if( (this->points[i]).first != NULL && (this->points[i].first)->app_evalue < 0 )
		{
			// Calling Pose constructor for the current chromosome
			Pose::Pose Pose((this->points[i]).first, i, this->order[i], this->reachDist[i], this->Population->Temperature, (this->points[i]).second);
			// OPTICS.push(Pose);
            this->OPTICS.push_back(Pose);
		}
	}
    std::sort(this->OPTICS.begin(), this->OPTICS.end(), PoseClassifier::PoseClassifier());

	// Build BindingModes (aggregation of Poses in BindingModes)
    int i = 0; // used to have an idea of the number of loop completed (iterators are less convenient for that information while debugging)
 	BindingMode::BindingMode Current(this->Population);
	for(std::vector< Pose >::iterator it = this->OPTICS.begin(); it != this->OPTICS.end(); ++i, ++it)
	{

        if(isUndefinedDist(it->reachDist) && it->order == 0) { Current.add_Pose(*it); }
        else if(it->reachDist <= this->FA->cluster_rmsd*(1 + RandomProjectedNeighborsAndDensities::sizeTolerance)) { Current.add_Pose(*it); }

       	if(it->reachDist > this->FA->cluster_rmsd*(1 + RandomProjectedNeighborsAndDensities::sizeTolerance) || isUndefinedDist(it->reachDist))
        {
			if(Current.get_BindingMode_size() >= this->minPoints)
			{
				this->Population->add_BindingMode(Current);
                Current.clear_Poses();
            }
        }
	}
}

void FastOPTICS::output_OPTICS(char* end_strfile, char* tmp_end_strfile)
{
	// output variables
    char sufix[25];

    // priting OPTICS variable to '__minPoints.optics' file
    sprintf(sufix,"__%d.optics",this->minPoints);
    strcpy(tmp_end_strfile,end_strfile);
	strcat(tmp_end_strfile,sufix);
    FILE* outfile;// = fopen(tmp_end_strfile,"w");
	if(!OpenFile_B(tmp_end_strfile,"w",&outfile))
	{
		Terminate(6);
	}
	else
	{
        fprintf(outfile, "#ORDER\t#INDEX\t#reachDist\t#CF\t#prevRMSD\t#nextRMSD\n");
		for(std::vector<Pose>::iterator it = this->OPTICS.begin(); it != this->OPTICS.end(); ++it)
		{
			float prevDist = (it == this->OPTICS.begin()) 	? 0.0f : this->compute_vect_distance(it->vPose, (it-1)->vPose);
			float nextDist = ((it+1) == this->OPTICS.end()) ? 0.0f : this->compute_vect_distance(it->vPose, (it+1)->vPose);
            if(!isUndefinedDist(it->reachDist))	fprintf(outfile, "%d\t%d\t%8g\t%8g\t%8g\t%8g\n", it->order, it->chrom_index, it->reachDist, it->CF, prevDist, nextDist);
            else fprintf(outfile, "%d\t%d\t%8g\t%8g\t%8g\t%8g\n", it->order, it->chrom_index, UNDEFINED_DIST, it->CF, prevDist, nextDist);
		}
   }
   CloseFile_B(&outfile,"w");;//fclose(outfile);
}

void FastOPTICS::output_3d_OPTICS_ordering(char* end_strfile, char* tmp_end_strfile)
{
	// File and Output variables declarations
    cfstr CF; /* complementarity function value */
    resid *pRes = NULL;
    cfstr* pCF = NULL;
	char sufix[25];
    char remark[MAX_REMARK];
	char tmpremark[MAX_REMARK];

	sprintf(sufix, "__%d.optics.pdb", this->minPoints);
	strcpy(tmp_end_strfile, end_strfile);
	strcat(tmp_end_strfile,sufix);
	FILE* outfile;
	if(!OpenFile_B(tmp_end_strfile, "w", &outfile))
	{
		Terminate(5);
	}
	else
	{
        int nModel = 1;
		for(std::vector<Pose>::iterator Pose = this->OPTICS.begin(); Pose != this->OPTICS.end(); ++Pose, ++nModel)
		{
			for(int k = 0; k < this->GB->num_genes; ++k) this->FA->opt_par[k] = Pose->chrom->genes[k].to_ic;
			CF = ic2cf(this->FA, this->VC, this->atoms, this->residue, this->cleftgrid, this->GB->num_genes, this->FA->opt_par);
			strcpy(remark,"REMARK optimized structure\n");
			sprintf(tmpremark,"REMARK Fast OPTICS clustering algorithm used to order Poses in OPTICS\n");
			strcat(remark,tmpremark);
			sprintf(tmpremark,"REMARK CF=%8.5f\n",get_cf_evalue(&CF));
			strcat(remark,tmpremark);
			sprintf(tmpremark,"REMARK CF.app=%8.5f\n",get_apparent_cf_evalue(&CF));
			strcat(remark,tmpremark);

			for(int j = 0; j < this->FA->num_optres; ++j)
			{
				pRes = &this->residue[this->FA->optres[j].rnum];
				pCF  = &this->FA->optres[j].cf;

		        sprintf(tmpremark,"REMARK optimizable residue %s %c %d\n", pRes->name, pRes->chn, pRes->number);
		        strcat(remark,tmpremark);

		        sprintf(tmpremark ,"REMARK CF.com=%8.5f\n", pCF->com);
		        strcat(remark, tmpremark);
		        sprintf(tmpremark ,"REMARK CF.sas=%8.5f\n", pCF->sas);
		        strcat(remark, tmpremark);
		        sprintf(tmpremark ,"REMARK CF.wal=%8.5f\n", pCF->wal);
		        strcat(remark, tmpremark);
		        sprintf(tmpremark ,"REMARK CF.con=%8.5f\n", pCF->con);
		        strcat(remark, tmpremark);
		        sprintf(tmpremark, "REMARK Residue has an overall SAS of %.3f\n", pCF->totsas);
		        strcat(remark, tmpremark);
			}

		    for(int j=0; j < this->FA->npar; ++j)
			{
				sprintf(tmpremark, "REMARK [%8.3f]\n",this->FA->opt_par[j]);
				strcat(remark,tmpremark);
			}

			// 4. if(REF) prints RMSD to REF
			if(this->FA->refstructure == 1)
			{
				bool Hungarian = false;
				sprintf(tmpremark,"REMARK %8.5f RMSD to ref. structure (no symmetry correction)\n",
				calc_rmsd(this->FA,this->atoms,this->residue,this->cleftgrid,this->FA->npar,this->FA->opt_par, Hungarian));
				strcat(remark,tmpremark);
				Hungarian = true;
				sprintf(tmpremark,"REMARK %8.5f RMSD to ref. structure     (symmetry corrected)\n",
				calc_rmsd(this->FA,this->atoms,this->residue,this->cleftgrid,this->FA->npar,this->FA->opt_par, Hungarian));
				strcat(remark,tmpremark);
			}

			// 5. write_pdb(FA,atoms,residue,tmp_end_strfile,remark)
			if(Pose == this->OPTICS.begin() && Pose+1 == this->OPTICS.end())
			{
				// case where there is only one pose to write (Pose == OPTICS.begin() && Pose++ == OPTICS.end())
				write_MODEL_pdb(true, true, nModel, this->FA,this->atoms,this->residue,tmp_end_strfile,remark);
			}
			else if(Pose == this->OPTICS.begin())
			{
				// first MODEL to be written
				write_MODEL_pdb(true, false, nModel, this->FA,this->atoms,this->residue,tmp_end_strfile,remark);
			}
			else if(Pose+1 == this->OPTICS.end())
			{
				// last MODEL to be written
				write_MODEL_pdb(false, true, nModel, this->FA,this->atoms,this->residue,tmp_end_strfile,remark);
			}
			else
			{
				// any MODEL in between to be written
				write_MODEL_pdb(false, false, nModel, this->FA,this->atoms,this->residue,tmp_end_strfile,remark);
			}
	}
}
}

std::vector<float> FastOPTICS::Vectorized_Chromosome(chromosome* chrom)
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
                    vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].coor[i] - this->FA->ori[i]);
					// vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].dis);
					// vChrom[i] *= vChrom[i];
				}
				if(i == 1)
				{
					vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].coor[i] - this->FA->ori[i]);
					// vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].ang);
					// vChrom[i] = static_cast<float>( RandomDouble( (*chrom).genes[j].to_int32) );
					// vChrom[i] *= vChrom[i];
				}
				if(i == 2)
				{
                    vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].coor[i] - this->FA->ori[i]);
					// vChrom[i] = static_cast<float>(this->cleftgrid[static_cast<unsigned int>((*chrom).genes[j].to_ic)].dih);
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

std::vector<float> FastOPTICS::Vectorized_Cartesian_Coordinates(int chrom_index)
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

void FastOPTICS::ExpandClusterOrder(int ipt)
{
    std::priority_queue< ClusterOrdering, std::vector<ClusterOrdering>, ClusterOrderingComparator::ClusterOrderingComparator > queue;
	ClusterOrdering tmp(ipt,0,1e6f);
	queue.push(tmp);

    while(!queue.empty())
	{
		ClusterOrdering current = queue.top();
		queue.pop();
		int currPt = current.objectID;

		if(this->processed[currPt] == true) continue;

		this->order[this->iOrder] = currPt;
		// incrementing STATIC rank ordering ()
		// FastOPTICS::iOrder++;
		this->iOrder++;
		this->processed[currPt] = true;

		float coredist = this->inverseDensities[currPt];
		for( std::vector<int>::iterator it = this->neighbors[currPt].begin(); it != this->neighbors[currPt].end(); ++it)
		{
			int iNeigh = *it;
			if(this->processed[iNeigh] == true) continue;

			float nrdist = this->compute_distance(this->points[iNeigh], this->points[currPt]);

			if(coredist > nrdist)
				nrdist = coredist;
			if(isUndefinedDist(this->reachDist[iNeigh]))
				this->reachDist[iNeigh] = nrdist;
			else if(nrdist < this->reachDist[iNeigh])
				this->reachDist[iNeigh] = nrdist;
			tmp = ClusterOrdering::ClusterOrdering(iNeigh, currPt, nrdist);
			queue.push(tmp);
		}
        std::make_heap(const_cast<ClusterOrdering*>(&queue.top()),
                       const_cast<ClusterOrdering*>(&queue.top() + queue.size()),
                       ClusterOrderingComparator::ClusterOrderingComparator()
                       );
	}
}

float FastOPTICS::compute_distance(std::pair< chromosome*,std::vector<float> > & a, std::pair< chromosome*,std::vector<float> > & b)
{
	float distance = 0.0f;

	// simple distance calculation below
	for(int i = 0; i < this->nDimensions; ++i)
	{
		float tempDist = (a.second[i]-b.second[i]);
        distance +=  tempDist * tempDist;
	}

   	return sqrtf(distance);
}

float FastOPTICS::compute_vect_distance(std::vector<float> a, std::vector<float> b)
{
	float distance = 0.0f;

	// simple distance calculation below
	for(int i = 0; i < this->nDimensions; ++i)
    {
		float tempDist = (a[i]-b[i]);
        distance +=  tempDist * tempDist;
    }

	return sqrtf(distance);
}

int FastOPTICS::get_minPoints() { return this->minPoints; }
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
        // std::vector<float> currentRp = this->Randomized_InternalCoord_Vector();
        std::vector<float> currentRp(this->Randomized_CartesianCoord_Vector());

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
        std::shuffle(projInd.begin(), projInd.end(), gen);

		int i = 0;
        for(std::vector<int>::iterator it = projInd.begin(); it != projInd.end(); ++it, i++)
		{
			int cind = (*it);
			// look this line to be sure that the vector is pushed in this->projectedPoints
			this->projectedPoints[cind] = tempProj[i];
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
		for(int l = 0; l < nElements-splitPos; ++l)
		{
			 ind3[l] = ind[l+splitPos];
//			ind3.push_back(ind[l+splitPos]);
		}
		this->SplitUpNoSort(ind3,dim+1);
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
				cneighs.push_back(oldind);
				std::sort(cneighs.begin(),cneighs.end());
				cneighs.erase(std::unique(cneighs.begin(), cneighs.end()), cneighs.end());
			}

			std::vector<int> & cneighs2 = neighs.at(oldind);
			// itPos = std::find(cneighs2.begin(), cneighs2.end(), ind);
			// if(itPos == cneighs2.end()) // element not found in cneighs2
			if( !std::binary_search(cneighs2.begin(), cneighs2.end(), ind) )
			{
				cneighs2.push_back(oldind);
				std::sort(cneighs2.begin(),cneighs2.end());
				cneighs2.erase(std::unique(cneighs2.begin(), cneighs2.end()), cneighs2.end());
			}
		}

	}
}

void FastOPTICS::normalizeDistances()
{
	float max = 0.0f;
	std::vector<float> & Distances = this->reachDist;
	for(std::vector<float>::iterator it = Distances.begin(); it != Distances.end(); ++it) if(*it > max && !isUndefinedDist(*it)) max = *it;
	for(std::vector<float>::iterator it = Distances.begin(); it != Distances.end(); ++it) if(!isUndefinedDist(*it)) *it /= max;
}

std::vector<float> RandomProjectedNeighborsAndDensities::Randomized_InternalCoord_Vector()
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
					vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].coor[i] - this->top->FA->ori[i]);
//					vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].dis);
//                    vChrom[i] *= vChrom[i];
				}
				if(i == 1)
				{
                    vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].coor[i] - this->top->FA->ori[i]);
//					vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].ang);
//                    vChrom[i] *= vChrom[i];
				}
				if(i == 2)
				{
                    vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].coor[i] - this->top->FA->ori[i]);
//					vChrom[i] = static_cast<float>(this->top->cleftgrid[static_cast<unsigned int>(doubleGeneIC)].dih);
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

std::vector<float> RandomProjectedNeighborsAndDensities::Randomized_CartesianCoord_Vector()
{
    float norm = 0.0f;

    int i = 0,j = 0,l = 0,m = 0;
	int cat;
	int rot;

	uint grd_idx;
	int normalmode=-1;
	int rot_idx=0;

    std::vector<float> vChrom(this->nDimensions);

	int npar = this->top->GB->num_genes;

	for(i=0;i<npar;i++){ this->top->FA->opt_par[i] = genetoic(&this->top->gene_lim[i],roll_die()); }

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

	if(normalmode > -1) alter_mode(this->top->atoms,this->top->residue,this->top->FA->normal_grid[normalmode],this->top->FA->res_cnt,this->top->FA->normal_modes);

	/* rebuild cartesian coordinates of optimized residues*/
    for(i=0;i<this->top->FA->nors;i++) buildcc(this->top->FA,this->top->atoms,this->top->FA->nmov[i],this->top->FA->mov[i]);

	// residue that is optimized geometrically (ligand)
	l=this->top->atoms[this->top->FA->map_par[0].atm].ofres;

	rot=this->top->residue[l].rot;
    m=0;
	for(i=this->top->residue[l].fatm[rot];i<=this->top->residue[l].latm[rot];i++)
	{
		for(j=0;j<3;j++)
		{
			vChrom[m*3+j] = this->top->atoms[i].coor[j];
			norm += vChrom[m*3+j]*vChrom[m*3+j];
		}
        ++m;
	}
	norm = sqrtf(norm);
//   	for(i = 0; i < this->nDimensions; ++i) vChrom[i] /= norm;
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

void RandomProjectedNeighborsAndDensities::output_projected_distance(char* end_strfile, char* tmp_end_strfile)
{
	char sufix[25];
	sprintf(sufix, "__%d.projDist", this->top->minPoints);
	strcpy(tmp_end_strfile, end_strfile);
	strcat(tmp_end_strfile,sufix);
	FILE* outfile;
	if(!OpenFile_B(tmp_end_strfile,"w",&outfile))
	{
		Terminate(5);
	}
	else
	{
		for(int i = 0; i < this->N; ++i)
		{
            for(int j = 0; j < this->nProject1D; ++j)
			{
                std::vector<float> & it = this->projectedPoints.at(j);
				fprintf(outfile, "%6f", it[i]);
                if(j < this->nProject1D-1) fprintf(outfile, "\t");
                else fprintf(outfile, "\n");
}
		}
	}
	CloseFile_B(&outfile,"w");;//fclose(outfile);
}

int roll_die()
{
    std::uniform_int_distribution<> dist(0, MAX_RANDOM_VALUE);
    return dist(gen);
}
int roll_rand_die()
{
	int n;
	return rand()%n;
}