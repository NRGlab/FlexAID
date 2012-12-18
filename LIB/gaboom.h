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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#include "flexaid.h"
#include "Vcontacts.h"

#define MAX_NUM_GENES 100
#define MAX_NUM_CHROM 1000
#define MAX_GEN_LENGTH 32                         // in number of bits
#define MAX_RANDOM_VALUE 2147483647              // upper bound of 32-bit integer
#define SAVE_CHROM_FRACTION 1.0

using namespace std;

typedef boost::mt19937 RNGType;

struct genelimits_struct{
	double max;
	double min;
	double del;
	double bin;
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
	double evalue;   // function evaluation of chromosome
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
	int	     print_int;

	char         pop_init_method[9];
	char         pop_init_file[20];
	char         fitness_model[9];
	char         rep_model[9];
	
};
typedef struct GB_Global_struct GB_Global;

/***********************************************************************/
/*        1         2         3         4         5         6          */
/*234567890123456789012345678901234567890123456789012345678901234567890*/
/*        1         2         3         4         5         6         7*/
/***********************************************************************/
int   GA(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome** chrom,chromosome** chrom_snapshot,genlim** gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid,char gainpfile[], int* memchrom, double (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int, double*));
int   check_state(char* pausefile, char* abortfile, char* stopfile, int interval);
void  crossover(gene *john,gene *mary,int num_genes);
void  mutate(gene *john,int num_genes,double mut_rate);
void  bin_print(int dec,int len);
void  read_gainputs(FA_Global* FA,GB_Global* GB,genlim** gene_lim,int*,int*,char file[]);

void   set_bins(genlim* gene_lim);
double set_bins(genlim* gene_lim, int num_genes);
double genetoic(const genlim* gene_lim, boost::int32_t gene);
int ictogene(const genlim* gene_lim, double ic);

int RandomInt(double frac);
double RandomDouble();
double RandomDouble(boost::int32_t dice);

void  swap_chrom(chromosome * x, chromosome * y);
void  quicksort_evalue(chromosome* list, int m, int n);
void  quicksort_fitnes(chromosome* list, int m, int n);
int   remove_dups(chromosome* list, int num_chrom, int num_genes);

void  populate_chromosomes(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim, atom* atoms,resid* residue,gridpoint* cleftgrid,char method[], double (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*), char file[], int offset, int print, boost::variate_generator< RNGType, boost::uniform_int<> > &);
double eval_chromosome(FA_Global* FA,GB_Global* GB,VC_Global* VC,const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,gene* john, double (*function)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*));
void  calculate_fitness(FA_Global* FA,GB_Global* GB,VC_Global* VC,chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,char method[],int pop_size, int print,double (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int, double*));
void  reproduce(FA_Global* FA,GB_Global* GB,VC_Global* VC, chromosome* chrom,const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,char rmodel[], double mutprob, double crossprob, int print, boost::variate_generator< RNGType, boost::uniform_int<> > &,double (*target)(FA_Global*,VC_Global*,atom*,resid*,gridpoint*,int,double*));
void  print_pop(const chromosome* chrom,const genlim* gene_lim,int numc, int numg);
void  print_chrom(const chromosome* chrom, int num_genes, int real_flag);
void  print_chrom(const gene* genes, int num_genes, int real_flag);
void  print_par(const chromosome* chrom,const genlim* gene_lim,int num_chrom,int num_genes);
int   roullete_wheel(const chromosome* chrom,int n);
int   cmp_chrom2pop(chromosome* chrom,chromosome* c, int num_genes,int start, int last);
int   cmp_chrom2pop_int(const chromosome* chrom,const gene* genes, int num_genes,int start, int last);
int   cmp_chrom2rotlist(psFlexDEE_Node psFlexDEE_INI_Node, const chromosome* chrom, const genlim* gene_lim,int gene_offset, int num_genes, int tot, int num_nodes);
int   cmp_chrom2pop(const chromosome* chrom,const gene* genes, int num_genes,int start, int last);
void  save_snapshot(chromosome* chrom_snapshot, const chromosome* chrom, int num_chrom, int num_genes);
void  cluster(FA_Global* FA, GB_Global* GB, VC_Global* VC,chromosome* chrom, genlim* gene_lim, atom* atoms, resid* residue,gridpoint* cleftgrid, int memchrom, char* end_strfile, char* tmp_end_strfile, char* dockinp, char* gainp);
//long long time_seed();
double calc_rmsp(int npar, const gene* g1, const gene* g2);
void write_par(const chromosome* chrom,const genlim* gene_lim,int ger, char* outfile,int num_chrom,int num_genes);
void adapt_prob(GB_Global* GB,double fitnes1,double fitnes2, double* mut_prob, double* cross_prob);
void fitness_stats(GB_Global* GB, const chromosome* chrom,int nchrom);
float  calc_rmsd_chrom(FA_Global* FA,GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,int npar, int chrom_a, int chrom_b); // calculates RMSD between chromossomes
int    write_rrg(FA_Global* FA,GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue, gridpoint* cleftgrid, char* outfile);        // writes GA output during simulation
int    write_rrd(FA_Global* FA,GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,int* Clus_GAPOP,float* Clus_RMSDT,char outfile[]);                            // writes output of population from GA
void   partition_grid(FA_Global* FA,chromosome* chrom,genlim* gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid,int pop_size, int expansion);        // partition grid size where favorable conformations are found
void   slice_grid(FA_Global* FA,genlim* gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid);                      // slice grid symmetrically in half


#endif // include guard
