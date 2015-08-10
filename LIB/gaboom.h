#ifndef GABOOM_H
#define GABOOM_H

#ifdef _WIN32
# ifndef _CRT_SECURE_NO_WARNINGS
#  define _CRT_SECURE_NO_WARNINGS
# endif
#endif

#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "boost/cstdint.hpp"
#include "boost/math/special_functions/fpclassify.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <limits.h>

#include "flexaid.h"
#include "Vcontacts.h"

#define MAX_NUM_GENES 100
#define MAX_NUM_CHROM 1000
#define MAX_GEN_LENGTH 32                         // in number of bits
#define MAX_RANDOM_VALUE 2147483647              // upper bound of 32-bit integer
#define SAVE_CHROM_FRACTION 1.0

#define QS_TYPE double
#define QS_ASC(a,b) ((a)-(b))
#define QS_DSC(a,b) ((b)-(a))
#define K(i,j,n) (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1
//#define K(i,j,n) ( (i < j) ? (i*n+j) : (j*n+i) )

using namespace std;

typedef boost::mt19937 RNGType;

struct genelimits_struct{
	double max;
	double min;
	double del;
	double bin;
	double nbin;
	int map;              // mapping gene (maps into an array)
};
typedef struct genelimits_struct genlim;

struct gene_struct{
	boost::int32_t  to_int32;
	double          to_ic;
};
typedef struct gene_struct gene;

struct chromosome_struct{
	gene*  genes;    // pointer to genes array
	cfstr  cf;       // function evaluation of chromosome
	double evalue;   // NOT the apparent cf
	double app_evalue;   // THE apparent cf
	double fitnes;   // fitness score of chromosome
	char   status;   /* status, n -> eval is correct 
			    o -> need to recalculate eval
			 */
};
typedef struct chromosome_struct chromosome;

struct GB_Global_struct{
	//long long    seed;

	int          num_chrom;
	int          num_genes;
	int          max_generations;
	
	double        alpha;
	double        peaks;
	double        scale;
	double        sig_share;

	double        mut_rate;
	double        cross_rate;
	double        ini_mut_prob;
	double        end_mut_prob;

	double        pbfrac;
	int           ssnum;
	int           intragenes;
	
	// adapt Genetic Algorithm operators (mutation,crossover)
	int          adaptive_ga;
	double       fit_avg;
	double       fit_max;
	double       fit_low;
	double       fit_high;
	double       k1,k2,k3,k4;

	int          outgen;
	int          rrg_skip;

	int          num_print;
	int	     	print_int;

	char         pop_init_method[9];
	char         pop_init_file[MAX_PATH__];
	char         fitness_model[9];
	char         rep_model[9];
	int          duplicates;
    
};
typedef struct GB_Global_struct GB_Global;

// Density Peaks clustering algorithm data structure definitions
struct ClusterChrom
{
	uint index;						// original index in chrommose* chrom (at the input of density_cluster function)
	bool isHalo;					// is part of cluster core?
	bool isCenter;					// is cluster center?
	bool isBorder;					// is part of the border region
	bool isClustered;
	chromosome* Chromosome;			// Chromosomes list
	int Cluster;					// Assigned Cluster
	int Density;					// Density of points in distance cut-off
	double CF;						// Complementarity Function value
	float PiDi;						// Density x Distance
	float Distance;					// Nearest highest density peak distance
	float Coord[3*MAX_ATM_HET];		// Cartesian Coordinates
	struct ClusterChrom* DP;		// Nearest Density Peak (point of higher density)
}; typedef struct ClusterChrom ClusterChrom;

struct Cluster_struct
{
       int ID;						 // assigned cluster number (ID)
       int Frequency;				 // observation frequency of this cluster (number of representatives in cluster)
       double totCF;
       double lowestCF;				
       // Pointer to best CF value in cluster
       ClusterChrom* Representative; // Pointer to the ClusterChrom individual with the lowest CF in cluster
       // Pointer to Cluster Center
       ClusterChrom* Center;		 // Queue of ClusterChrom (first element is the cluster center)
};
typedef struct Cluster_struct DPcluster;
/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int   GA(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome** chrom,chromosome** chrom_snapshot,genlim** gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid,char gainpfile[], int* memchrom, cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int, double*));
int   check_state(char* pausefile, char* abortfile, char* stopfile, int interval);
void  QuickSort(chromosome*, int, int, bool);
void  QuickSort_Clusters(int*, int*, double*, double*, int*, int, int);
void  swap_clusters(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*);
void  crossover(gene *john,gene *mary,int num_genes, int intragenes);
void  mutate(gene *john,int num_genes,double mut_rate);
void  bin_print(int dec,int len);
void  read_gainputs(FA_Global* FA,GB_Global* GB,int*,int*,char file[]);
int   deelig_search(struct deelig_node_struct* root_node, int* deelig_list, int fdih);
int   filter_deelig(FA_Global* FA, GB_Global* GB, chromosome* chrom, gene* genes, int ci, atom* atoms, const genlim* gene_lim,
		   boost::variate_generator< RNGType, boost::uniform_int<> > & dice);

void   set_gene_lim(FA_Global* FA, GB_Global* GB, genlim* gene_lim);
long int read_pop_init_file(FA_Global* FA, GB_Global* GB, genlim* gene_lim, char* pop_init_file);
void  set_bins(genlim* gene_lim);
void  set_bins(genlim* gene_lim, int num_genes);
double calc_poss(genlim* gene_lim, int num_genes);
void validate_dups(GB_Global* GB, genlim* gene_lim, int num_genes);
double genetoic(const genlim* gene_lim, boost::int32_t gene);
int ictogene(const genlim* gene_lim, double ic);

int 	RandomInt(double frac);
double 	RandomDouble();
double 	RandomDouble(boost::int32_t dice);

void  	swap_chrom(chromosome * x, chromosome * y);
void  	copy_chrom(chromosome* dest, const chromosome* src, int num_genes);
int   	remove_dups(chromosome* list, int num_chrom, int num_genes);

FILE* 	get_update_file_ptr(FA_Global* FA);
void 	close_update_file_ptr(FA_Global* FA, FILE* outfile_ptr);
string 	generate_sig(gene genes[], int num_genes);

void  	generate_random_individual(FA_Global* FA, GB_Global* GB, atom* atoms, gene* genes, const genlim* gene_lim,
				 boost::variate_generator< RNGType, boost::uniform_int<> > &, int from_gene, int to_gene);
void  	populate_chromosomes(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim, atom* atoms,resid* residue,gridpoint* cleftgrid,char method[], cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*), char file[], long int at, int offset, int print, boost::variate_generator< RNGType, boost::uniform_int<> > &, map<string, int> &);
cfstr 	eval_chromosome(FA_Global* FA,GB_Global* GB,VC_Global* VC,const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,gene* john, cfstr (*function)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*));
void  	calculate_fitness(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,char method[],int pop_size, int print, cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int, double*));
int   	reproduce(FA_Global* FA,GB_Global* GB,VC_Global* VC, chromosome* chrom,const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,char rmodel[], double mutprob, double crossprob, int print, boost::variate_generator< RNGType, boost::uniform_int<> > &,map<string, int> &, cfstr (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*));
void  	print_pop(const chromosome* chrom,const genlim* gene_lim,int numc, int numg);
void  	print_chrom(const chromosome* chrom, int num_genes, int real_flag);
void  	print_chrom(const gene* genes, int num_genes, int real_flag);
void  	print_par(const chromosome* chrom,const genlim* gene_lim,int num_chrom,int num_genes, FILE* outfile_ptr);
int   	roullete_wheel(const chromosome* chrom,int n);
int   	cmp_chrom2pop(chromosome* chrom,chromosome* c, int num_genes,int start, int last);
int   	cmp_chrom2pop_int(const chromosome* chrom,const gene* genes, int num_genes,int start, int last);
int   	cmp_chrom2rotlist(psFlexDEE_Node psFlexDEE_INI_Node, const chromosome* chrom, const genlim* gene_lim,int gene_offset, int num_genes, int tot, int num_nodes);
int   	cmp_chrom2pop(const chromosome* chrom,const gene* genes, int num_genes,int start, int last);
void  	save_snapshot(chromosome* chrom_snapshot, const chromosome* chrom, int num_chrom, int num_genes);
void  	cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC,chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue,gridpoint* cleftgrid, int memchrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp);
void  	DensityPeak_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gen_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int memchrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp);
void 	FastOPTICS_cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC, chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, int nChrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp);
//long long time_seed();
double 	calc_rmsp(int npar, const gene* g1, const gene* g2, const optmap* map_par, gridpoint* cleftgrid);
void 	write_par(const chromosome* chrom,const genlim* gene_lim,int ger, char* outfile,int num_chrom,int num_genes);
void 	adapt_prob(GB_Global* GB,double fitnes1,double fitnes2, double* mut_prob, double* cross_prob);
void 	fitness_stats(GB_Global* GB, const chromosome* chrom,int nchrom);
float  	calc_rmsd_chrom(FA_Global* FA,GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,int npar, int chrom_a, int chrom_b, float*, float*, bool calc_rmsd); // calculates RMSD between chromossomes
int    	write_rrg(FA_Global* FA,GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue, gridpoint* cleftgrid, char* outfile);        // writes GA output during simulation
int    	write_rrd(FA_Global* FA,GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,int* Clus_GAPOP,float* Clus_RMSDT,char outfile[]);   
int 	write_DensityPeak_rrd(FA_Global* FA, GB_Global* GB, const chromosome* chrom, const genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, ClusterChrom* Chrom, DPcluster* Clust, float* RMSD, char outfile[]);
void   	partition_grid(FA_Global* FA,chromosome* chrom,genlim* gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid,int pop_size, int expansion);        // partition grid size where favorable conformations are found
void   	slice_grid(FA_Global* FA,genlim* gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid);                      // slice grid symmetrically in half

// Density Peaks Clustering algorithm function declarations
void 	QuickSort_Cluster_by_CF(DPcluster* Clust, bool Entropic, int beg, int end);
void 	swap_clusters(DPcluster* xClust, DPcluster* yClust);
float 	getDistanceCutoff(float* RMSD, int num_chrom);
void 	QuickSort_ChromCluster_by_CF(ClusterChrom* Chrom, int num_chrom, int beg, int end);
void 	QuickSort_ChromCluster_by_higher_Density(ClusterChrom* Chrom, int num_chrom, int beg, int end);
void 	QuickSort_ChromCluster_by_lower_Density(ClusterChrom* Chrom, int num_chrom, int beg, int end);
void 	swap_elements(ClusterChrom* Chrom, ClusterChrom* ChromX, ClusterChrom* ChromY, int num_chrom);
float 	calculate_stddev(ClusterChrom* Chrom, int num_chrom);
float 	calculate_mean(ClusterChrom* Chrom, int num_chrom);
int 	DistanceComparator(const void*, const void*);

#endif // include guard
