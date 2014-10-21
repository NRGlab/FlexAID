// FLEXAID HEADER FILE

#ifndef FLEXAID_H
#define FLEXAID_H

#ifdef _WIN32

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <utility>
#include <map>

//#include "/home/francis/valgrind/include/valgrind/memcheck.h"

//const int endian_t = 1;
//#define IS_BIG_ENDIAN() ( ( *(char *) &endian_t ) == 0 ) // cross-platform development

#define MAX_ENERGY_POINTS 5000      // probability function distribution
#define MAX_CONTRIBUTIONS 10
#define MAX_SHORTEST_PATH 25        // max number of atom to reach any atom of the same molecule
#define MAX_PATH__ 250              // max size of path length
#define MAX_REMARK 5000             // max size of comment length
#define MAX_NUM_RES 250             // max number of residues allowed 
#define MAX_NUM_MODES 5             // max number of normal modes
#define MAX_NUM_ATM 10000           // max number of atoms allowed 
#define MAX_ATM_HET 200             // max number of atoms in het groups allowed 
//#define MAX_FLEX_BONDS 20           // max number of flexible bonds for het groups
#define MAX_GRID_POINTS 1000        // Total number of points to anchor ligand
#define MAX_NORMAL_GRID 5000        // max number of grid points in normal grid
#define MAX_SPHERE_POINTS 610       // Total number of points in atom surface sphere 
#define MAX_OPT_RES 1               // max number of residues to be optimized 
#define POLYHEDRON_PENALTY 1e6f     // skipped individuals WAL term penaly
#define CLASH_THRESHOLD 1e4         // individuals that do not pass the clash filter of Vcontacts
#define MBNDS 4                     // max number of cov bonds an atom can have 
#define MAX_BONDED 20               // max number of bonded atom list
#define MAX_PAR 100                 // max number of parameters for simplex optimization, not used at the moment
#define MAX_ROTLIBSIZE 500
#define MAXFLXSC 100
#define MAX_CLOSE_DIST 10
#define RMSD_THRESHOLD 2.0

#define KWALL    1.0e6
#define KANGLE   1.0e2
#define KDIST    1.0e3
#define DEE_WALL_THRESHOLD 50.0

#define PI 3.141592
#define E  2.718281
#define Rw 1.4f

#define NEW(p,type)     if ((p=(type *) malloc (sizeof(type))) == NULL) { \
		printf ("Out of Memory!\n");				\
		exit(0);						\
	}

#define FREE(p)         if (p) { free ((char *) p); p = NULL; }


#define ADD( head, p )  if ( head )  {		\
		p->next = head;			\
		p->prev = head->prev;		\
		head->prev = p;			\
		p->prev->next = p;		\
	}					\
	else {					\
		head = p;			\
		head->next = head->prev = p;	\
	}

/*
  #define DELETE( head, p ) if ( head )  { \
  if ( head == head->next ) \
  head = NULL;  \
  else if ( p == head ) \
  head = head->next; \
  s                                p->next->prev = p->prev;  \
  p->prev->next = p->next;  \
  FREE( p ); \
  }
*/

typedef unsigned int uint;

struct cf_str{  // Complementarity Function value structure
	double com;    // complementarity value
	double con;    // constraint value
	double wal;    // wall term
	double sas;    // solvent accessibility surface
	double totsas; // overall sas of molecule
	int   rclash; // flag that shows whether the residue is making steric clashes
};
typedef struct cf_str cfstr;

struct grid_point_struct{ //fixed points in space internal coordinates
	int   index;   // index needed to sort in order of distance
	int   number;  // grid point C atom number
	float coor[3]; // coordinates of one intersection of the grid
	float dis;     // distance to PCG + (1,0,0)
	float ang;     // angle to PCG
	float dih;     // dihedral angle to PCG
};
typedef struct grid_point_struct gridpoint;

struct OptRes_struct{
	int     rnum;      // internal residue number
	int     type;      // type 0: protein - 1: ligand
	int     tot;       // number of atoms in res
	cfstr   cf;        // cf value
};
typedef struct OptRes_struct OptRes;

struct energy_values {
	float x;
	float y;
	struct energy_values* next_value;
};

struct energy_matrix {
	int type1;
	int type2;
	int weight;        // weights surface in contact vs. probability functions
	struct energy_values* energy_values;
};

struct constraint_str{
	int  rnum1;
	int  rnum2;
	char rnam1[5];
	char rnam2[5];
	char chn1;
	char chn2;
	int  anum1;
	int  anum2;
	int  inum1;
	int  inum2;
    
	float max_ang;    // threshold critical angle value
	float max_dist;   // threshold critical distance value

	int  id;          // constraint id
	int  type;        // type of constraint 0: interaction - 1:covalent bond
	float bond_len;   // minimal distance between 2 atoms before WAL term is calculated

	float interaction_factor;
	int   force_interaction;
};
typedef struct constraint_str constraint;

struct optmap_struct{  // optimization residues structure
	int typ; // type of atom to be optimized
	int atm; // number of atom to be optimized
	int bnd; // which flexible bond needs to be optimized
};
typedef struct optmap_struct optmap;

struct atom_struct{  // atom structure
	float  coor[3];     // processing coordinates
	float* coor_ref;    // reference coordinates
	float  coor_ori[3]; // original coordinates
	
	int    number;  // atom number according to PDB file
	float  radius;  // atomic radius
	int    type;    // atom type
	int    bond[5]; // atom number (not according to PDB) of covalently bonded atoms, bond[0] gives the total
	int    ofres;   // residue number to which an atom belongs.
	char   recs;    // flag which labels an atom as flexible or rigid
	float  dis;     // distance from rec[0] atom
	float  ang;     // angle between atom, rec[0] and rec[1]
	float  dih;     // dihedral between atom, rec[0], rec[1] and rec[2]
	float  shift;   // gives the angle shift from that of rec[3]'s atom
	float  acs;     // accessible contact surface
	int    ncons;   // number of constraint for atoms
	int    isbb;    // atom is a backbone atom
	int    graph;   // id of graph atom belongs to (ligands only)
	
	optmap* par;    // if this atom defines a variable (translational/rotational or dihedrals)
	constraint** cons; // points to constraint , if NULL no constraint to atom
	OptRes* optres;  // pointer to optimised residue list
	float** eigen;   // eigen vectors
	
	int    rec[4];  // atom number to be used when reconstructing the atom coordinates from internal coordinates
	char   name[5]; // atom name
	char   element[3]; // element name
};
typedef struct atom_struct atom;

struct residue_struct{   // residue structure
	char  name[4];  // residue name
	char  chn;      // chain name
	char  ins;      // insertion
	int   number;   // residue number according to PDB
	int*  fatm;     // number of first atom for a given rotamer of given residue
	int*  latm;     // number of last atom for a given rotamer of given residue
	int** bonded;   // bonded list of atoms (i x j matrix)
	char*** shortpath;  // shortest path of atoms (i x j matrix)
	int*** shortflex; // list of flexible bonds between 2 atoms following the shortest path
	int   type;     // labels residues as protein or ligand residues
	int   fdih;     // number of flexible bonds
	int*  bond;     // atoms whose rec[3] defines each flexible bond
	int*  gpa;      // global positioning atoms
	int   ter;      // if residue is a terminating residue
	int   trot;     // total number of rotamers
	int   rot;      // specific rotamer number, 0 is PDB
};
typedef struct residue_struct resid;

struct flexible_sc_struct{ // flexible side chain structure
	char  name[4];  // residue name
	char  chn;      // chain name
	int   num;      // residue PDB number
	int   inum;     // internal residue number
	float prob;     // prob of changing conformation

	int   cflag;    // used to skip side-chains when all neighbours are rigid
	int*  close;    // list of internal residue numbers that will need to be tested for clash
	int   nclose;   // number of element in list
};
typedef struct flexible_sc_struct flxsc;

/*
  class DEELig_Node {
  public:
  void DEELig_Node(void);
  private:
  struct deelig_node_struct* parent;
  std::map<int, struct deelig_node_struct*> childs;
  };
*/

struct deelig_node_struct{
	struct deelig_node_struct* parent;
	std::map<int, struct deelig_node_struct*> childs;
};
//typedef struct deelig_node_struct deelig_node;

struct RotLib_struct{ // Rotamer library entry records
	char  res[4];  // residue name
	int   nid;     // numeric id of rotamer
	char  name[9]; // rotamer name
	int   tot;     // total number of residue of that type
	int   obs;     // number of observations
	float pro;     // rotamer probability
	int   nchi;    // number of dihedral angles
	float chi[4];  // values of dihedral angles
	int   numrng;  // number of ranges
	float lowr[4][2]; // low-bound values in range
	float hghr[4][2]; // high-bound values in range
	float gaus[4][2]; // positive~negative 1/2Width 1/2Height for Gaussian-Poisson rotamer dist
};
typedef struct RotLib_struct rot;


typedef struct sphere_struct sphere;
struct sphere_struct {
	float radius;
	float center[3];
	sphere* prev;
};


typedef struct FlexDEE_Node_struct sFlexDEE_Node;
typedef sFlexDEE_Node* psFlexDEE_Node;

struct FlexDEE_Node_struct {
	int* rotlist;
  
	psFlexDEE_Node first;
	psFlexDEE_Node last;

	psFlexDEE_Node next;
	psFlexDEE_Node prev;
};


struct FA_Global_struct{
	optmap* map_par;                 // array of structure of mapping of optimization parameters

	int map_par_flexbond_first_index; // index in array of first flexbond map par
	optmap* map_par_flexbond_first;  // first map representing the flexible bonds
	optmap* map_par_flexbond_last;   // last map representing the flexible bonds

	int map_par_sidechain_first_index; // index in array of first flexbond map par
	optmap* map_par_sidechain_first; // first map representing the side-chain rotamers
	optmap* map_par_sidechain_last;  // last map representing the side-chain rotamers

 	double*  opt_par;                // optimization parameters  
	struct deelig_node_struct* deelig_root_node;   // termination criteria for ligand flexible
	int   deelig_flex;

	int    npar;                         // number of parameters
	int*  num_atm;                       // PDB num --> internal num mapping
	float* contributions;                // contributions of interactions from the complexe
	
	//atom  ori_ligatm[100];               // array to carry original atomic coordinates of ligand
	//int   num_ligatm;                    // number of atoms in original ligand
    
	//atom  atoms[MAX_NUM_ATM];            // array of atom structure
	int   atm_cnt;                       // total number of atoms including all rotamers
	int   atm_cnt_real;                  // total number of atoms real
	int   nflexbonds;                    // number of ligand flexible bonds
	
	int   normalize_area;                // normalize contact areas as a function of total surface area
	int   useacs;                        // normalize interactions by accessible contact surface
	float acsweight;                     // weighting factor
	
	int   is_protein;                    // PDBNAM is a protein molecule
	resid* resligand;                    // pointer to the ligand
	//int   is_nucleicacid;              // PDBNAM is a DNA/RNA molecule

	float ori[3];                        // coordinates of center of geometry of protein (PCG)
	float maxdst;                        // max distance from a protein atom to PCG 
	float cluster_rmsd;                  // rmsd between poses when clustering

	float permeability;                  // allow permeability or not between atoms
	float rotamer_permeability;          // rotamer acceptance vdw permeability
	int   intramolecular;                // consider intramolecular forces (ligand only)
	float solventterm;                   // solvent penalty term
	float intrafraction;                 // intramolecular fraction interaction

	constraint* constraints;             // list of constraints
	int num_constraints;                 // constraints counter
	float interaction_factor;            // interaction constraint factor
	int force_interaction;               // forces interaction although complementarity is negative
  
	int   max_results;                   // maximum result file(s) generated after GA.

	//resid residue[MAX_NUM_RES];          // array of residues structure
	int   res_cnt;                       // total number of residues
	int   exclude_het;                   // exclude HET groups when calculating CF
	int   remove_water;                  // exclude water molecules (HET=HOH)
	int   output_range;                  // outputs Sphere or Grid file(s)

	int     bloops;                      // exclude interactions with atoms n bloops away (exclude dis-ang preferably)
	OptRes* optres;                      // residues optimised (side-chains + ligand) for which CF needs to be calculated
	int     num_optres;                  // total number of residues optimised

	int   rotobs;                        // use rotamer observations, otherwise default Lovell's LIBrary
	int   rotout;                        // output rotamers in rotamers.pdb as pdb models
	
	int   num_het;                       // number of hetero groups read
	int   num_het_atm;                   // number of hetero atoms
	int   num_grd;                       // number of vertices in the grid of the cleft
	int   het_res[2];                    // table of residue number for each HET group 
	int   opt_grid;                      // flag determining whether to partition/splice the grid during GA
  
	int   nmov[2];                       // the next three items are used to determine
	int*  mov[2];                        // in which order to rebuilt atoms that are
	int   nors;                          // flexible
	int   opt_res[2];                    // list of residue numbers being optimized
  
	float spacer_length;                 // space length between intersections of the grid
  
	//atom* atoms_rmsd;                    // atoms to calculate rmsd
	//resid res_rmsd[2];                   // residues to calculate rmsd
	//int   natoms_rmsd;                   // number of atoms in atoms_rmsd
	//int   hrnum;

	int output_scored_only;              // ouptuts the ligand coordinates only in the results file
	int score_ligand_only;              // scores the ligand only despite sidechains are enabled.

	char vcontacts_self_consistency[6];  // A --> B and B --> A contacts self consistency
	char vcontacts_planedef;             // plane definition for vcontacts
    
	int   nrg_suite;                     // flag indicating if nrg_suite is enabled
	int   nrg_suite_timeout;             // specifies the maximum time for the suite to update the visuals (in seconds)
	int   translational;                 // flag indicating if translation degrees of freedom are enabled
	int   refstructure;                  // reference structure for rmsd calculation
  
	int* contacts;                       // matrix used for not calculating the same interaction twice
	struct energy_matrix* energy_matrix;        // potential energy parameters
	int   ntypes;	                     // number of atom types
	int   tspoints;                      // actual number of sphere points
	//long long seed_ini;                // initial seed for random num generator
	int   recalci;                       // recalculations counter (VCT)
	int   skipped;                       // atoms skipped due to faliures in generating the polyhedron
	int   clashed;                       // skipped individuals due to steric clashes
	int   omit_buried;                   // skip buried atoms in the Vcontacts procedure
	int   vindex;                        // use indexed boxes and atoms in Vcontacts index_proteins

	//rot    rotamer[MAX_ROTLIBSIZE];       // array of rotamer library rotamers OR observed rotamer list
	int    rotlibsize;                    // number of rotamers
	flxsc* flex_res;                      // list of flexible residues
	int    nflxsc;                        // number of flexible residues defined and found within the protein
	int    nflxsc_real;                   // number of flex. resi that have at least one accepted rotameric conformation
	int    pbloops;                       // Number of different binding site conformations generated

	int    useflexdee;                    // use dead-end elimination for flexible side-chains
	float  dee_clash;
	psFlexDEE_Node psFlexDEENode;         // starting Node in DEAD_END_ELIMINATION for side-chains  
	int    FlexDEE_Nodes;                 // number of Nodes  

	//int    normal_mode;                        // flag that enables normal mode
	int      normal_modes;                       // Number of Normal Modes combined   
	int      supernode;                          // flag for activation supernode mode
	int      normal_nodes;                       // Number of Nodes (C,CA,N atoms)
	int      normal_grid_points;                 // Number of grid points in normal grid

	float**  normal_grid;                        // 2-dimensional grid containing amplitudes of eigenvectors
	float**  eigenvector;                        // 2-dimensional grid containing eigen vectors
	int      num_eigen;
	
	double normalindex_min;                // boundaries of IC for normal mode in GA
	double normalindex_max;                // ...

	double*  del_opt_par;         // delta optimization parameters
	double*  min_opt_par;         // min value optimization parameters
	double*  max_opt_par;         // max value optimization parameters
	int*     map_opt_par;         // parameter needs to be mapped to an array

	double delta_angstron;               // delta in angstrons for optimization or move of distances
	double delta_angle;                  // delta in degrees for optimization or move of bond angle
	double delta_dihedral;               // delta in degrees for optimization or move of dihedral angle
	double delta_flexible;               // delta in degrees for optimization or move of flexible dihedrals (ligands)
	double delta_index;                  // delta in integer for optimization of indexes for spheres

	float sphere[MAX_SPHERE_POINTS][3];  // coordinates of the unit sphere

	char  metopt[3];                     // string for defining optimization method GA or other.
	char  bpkenm[3];                     // string for defining binding pocket enumeration (XS or PB).
	char  complf[4];                     // string for defining complementarity function (SPH or VCT)
	char  rngopt[7];                     // range of search
	resid rngres;                        // residue structure when search is around a given residue
	float rngcen[3];                     // coordinates of point around search is made
	float rngrad;                        // radius of search around point/residue
	float globalmin[3];                  // minimum coordinate value for protein atoms
	float globalmax[3];                  // maximum coordinate value for protein atoms
	float maxwidth;                      // maxwidth from min and max global coor of all protein atoms
	double dis_min,ang_min,dih_min;       // min values of distance, angle and dihedral
	double dis_max,ang_max,dih_max;       // max values of distance, angle and dihedral
	double index_min,index_max;           // max values of indexes of fixed spheres into space

	char  stp[2];                        // string for pausing output, NULL 

	char  rrgfile[MAX_PATH__];             // filename for reading initial pop for GA.

	char  base_path[MAX_PATH__];           // path leading to the executable of FlexAID
	char  dependencies_path[MAX_PATH__];   // path leading to dependencies files (e.g. AMINO.def, radii.dat, scr_bin.dat, ...)
	// if dependencies path is not set, dependencies will be searched in base_path unless forced

	char  state_path[MAX_PATH__];          // path leading to files .pause, .stop, .abort
	char  temp_path[MAX_PATH__];           // path where target.pdb is written, range files (grid/sphere), defaults to working dir

	// minimums used for dynamic allocation
	int    MIN_NUM_ATOM;
	int    MIN_NUM_RESIDUE;
	int    MIN_ROTAMER_LIBRARY_SIZE;
	int    MIN_ROTAMER;
	int    MIN_FLEX_BONDS;
	int    MIN_CLEFTGRID_POINTS;
	int    MIN_PAR;
	int    MIN_OPTRES;
	int    MIN_FLEX_RESIDUE;
	int    MIN_NORMAL_GRID_POINTS;
	int    MIN_CONSTRAINTS;
};
typedef struct FA_Global_struct FA_Global;


//-----------------------------------------------------------------------------------------
// PROTOTYPES:
//-----------------------------------------------------------------------------------------
float  assign_radius(char atm[]);                                  // assigns atomic radii
void   assign_radii(atom* atoms,resid* residue,int atmcnt);      // overrides default radii from FlexAID
void   assign_radii_types(FA_Global* FA, atom* atoms, resid* residue); // sets radius using the SYBYL types

// NEW GEOMETRY INLINE FUNCTIONS
float  dihedral(const float *,const float *,const float *,const float *);
float *cross_prod(float *, const float *, const float *);
float  dot_prod(const float *, const float *);
void   vec_sub(float *, const float *, const float *);
float  angle(const float *, const float *, const float *);
float  distance(const float *, const float *);
float  distance2(const float *, const float *);
float  zero(float, float, float);
float  dihang(float a[],float b[], float c[], float d[]);  // calculates dihedral angles
float  bndang(float a[],float b[], float c[]);             // calculates bond angles
float  sqrdist(float a[], float b[]);                      // calculates square distance
float  dist(float a[], float b[]);                     // calculates distance
float  distance_n(float a[], float b[], int n);            // calculates distance, n-dimensional
int    spfunction(FA_Global* FA,atom* atoms,resid*);                                // CF-SPHERE function

void   read_coor(FA_Global* FA,atom** atoms,resid** residue,char line[], char res_numold[]);          // reads PDB coordinates
void   read_conect(FA_Global* FA,atom** atoms,char line[]);                           // reads CONECT field for ligand
void   read_pdb(FA_Global* FA,atom** atoms,resid** residue, char* pdb_name);                          // reads PDB file and acts accordingly
void   read_grid(FA_Global* FA,gridpoint** cleftgrid,char file[]);                             // reads cleft sphere file
void   read_normalgrid(FA_Global* FA,char file[]);                       // reads normal mode grid file
void   read_eigen(FA_Global* FA,char file[]);                            // reads normal mode grid file
void   assign_eigen(FA_Global* FA,atom* atoms,resid* residue,int rescnt, int nmodes);  // assigns vector to atoms structure
void   alter_mode(atom *atoms, resid* residue, float *grid,int res_cnt, int nmodes);  // alter the backbone of the protein using normal mode

void   wif083(FA_Global* FA);                                            // generates surface points for cffunction
void   residue_conect(FA_Global* FA,atom* atoms,resid* residue,char file[]);                    // assigns covalent bonds for protein atoms
void   assign_types(FA_Global* FA,atom* atoms,resid* residue,char file[]);                      // assigns atom types for protein atoms
void   buildprob();                                        // build the rotamer probability list used as a roulette wheel
void   build_rotamers(FA_Global* FA,atom** atoms,resid* residue,rot* rotamer);         // build rotamer atoms in atoms structure
void   calc_cleftic(FA_Global* FA,gridpoint* cleftgrid);                                      // calculates dis, ang and dih for each dot(sphere) of the binding site
void   buildlist(FA_Global* FA,atom* atoms,resid* residue,int rnum, int bnum, int *tot, int lout[]);// creates list of atoms that need to be rebuilt
void   buildcc(FA_Global* FA,atom* atoms,int tot,int list[]);                        // creates cartesian coordinates from internal coords.
void   buildic(FA_Global* FA,atom* atoms,resid* residue,int rnum);                   // creates internal coordinates from cartesian.
void   add2_optimiz_vec(FA_Global* FA,atom* atoms,resid* residue,int val[], char chain, const char* extras);            // adds atoms that need to be optimized
void   realloc_par(FA_Global* FA, int* MIN_PAR); // reallocs memory for par in add2 function
void   read_lig(FA_Global* FA,atom** atoms,resid** residue,char ligfile[]);                           // reads ligands
void   read_input(FA_Global* FA,atom** atoms,resid** residue,rot** rotamer,gridpoint** cleftgrid,char input_file[]);                      // reads input file
void   print_surfmat(FA_Global* FA,float matrix[9][9], char* filename);
void   create_rebuild_list(FA_Global* FA,atom* atoms,resid* residue);                               // calls buildlist
void   bondedlist(atom* atoms,int anum, int nloops, int* nlist_ptr, int* blist, int* nnbr); // determines bonded atoms for cffunction
void   assign_shift(atom* atoms,resid* residue,int rnum, int natm, int *list, int **fdihlist);      // assigns shift atoms internally
int    write_pdb(FA_Global* FA,atom *atoms,resid* residue,char outfile[], char remark[]);           // writes PDB files
void   write_contributions(FA_Global* FA, FILE* outfile_ptr, bool positive);
const char* get_element(int type);
void   write_grid(FA_Global* FA, const gridpoint* cleftgrid,char filename[]);                        // writes starting grid to PDB viewable format with 'grid' ext. PYTHON function
void   calc_center(FA_Global* FA,atom* atoms,resid* residue);            // calculates center of geometry of protein
float  calc_rmsd(FA_Global* FA,atom* atoms,resid* residue, gridpoint* cleftgrid, int npar, const double* icv);       // calculates rmsd
void   ic_bounds(FA_Global* FA, char* rngopt);              // determines internal coordinates bounds for GA opt.
void   buildic_point(FA_Global* FA,float coor[],float *dis, float *ang, float *dih); // builds IC for one atom
void   buildcc_point(FA_Global* FA,float coor[], float dis, float ang, float dih);   // builds CC for one atom
int    read_rmsdst(FA_Global* FA,atom* atoms,resid* residue,char rmsdst_file[]);       // reads ref. structure input
void   write_sphere(FA_Global* FA,char filename[]);                                    // build sphere file 
void   read_rotlib(FA_Global* FA,rot** rotamer,char* libfile);                         // reads rotamer library
void   read_rotobs(FA_Global* FA,rot** rotamer,char* rotobsfile);                         // reads rotamer observation file
void   read_constraints(FA_Global* FA,atom* atoms, resid* residue, char* constfile);                         // reads constraints file
void   read_emat(FA_Global* FA, char* scr_bin_file);                 // reads interaction matrix
sphere* read_spheres(char filename[]);                               // reads the list of spheres from a cleft
gridpoint* generate_grid(FA_Global* FA, sphere* spheres, atom* atoms, resid* residue);            // builds the grid from loccen/locclf
void   assign_constraint_threshold(FA_Global* FA,atom* atoms,constraint* cons,int ncons);      // pre-calculate critical values

int    assign_constraint(FA_Global* FA,atom* atoms, resid* residue,constraint* cons);                         // assign constraint to atoms structure
void   update_constraint(atom* atoms, int index, constraint* cons);                                            // add new constraint in atoms structure
void   shortest_path(resid* residue, int tot, atom* atoms);
void   assign_shortflex(resid* residue, int tot, int fdih, atom* atoms);
void   update_bonded(resid* residue, int tot, int nlist, int* list, int* nbr);   // update bonded matrix for residue
void   update_optres(atom* atoms, resid* residue, int atm_cnt, OptRes* optres_ptr,int num_optres); // update atoms structure with optres pointers
void   set_intprob(flxsc* flex_res);                                         // sets internal probability for rotamer change
int    number_of_dihedrals(char res[]);                      // determines number of dihedral bonds for residues
void   read_flexscfile(FA_Global* FA,resid* residue,rot** rotamer,char* file, char* rotlib, char* rotobs);                        // reads flexible residues list
int    check_clash(FA_Global* FA,atom* atoms,resid* residue,int res_cnt,int total, int list[]);   // checks if there are clashes with rigid residues
void   build_close(FA_Global* FA, resid** residue, atom** atoms);

int  dee_pivot(psFlexDEE_Node psFlexDEEInsertNode,psFlexDEE_Node* psFlexDEENode, int sta, int end, int pivot, int pos, int num_rot);         // recursives through DEE list to add item in sorted way
void dee_first(psFlexDEE_Node Node, psFlexDEE_Node FirstNode);       // sets first record for all items
void dee_last(psFlexDEE_Node Node, psFlexDEE_Node LastNode);         // sets last record for all items
void dee_print(psFlexDEE_Node psFlexDEENode, int num_rot);           // prints DEE sorted list
int  cmp_FlexDEE_Nodes(psFlexDEE_Node Node1, psFlexDEE_Node Node2, int num_rot); // compares the DEE list between items

double get_apparent_cf_evalue(cfstr* cf);
double get_cf_evalue(cfstr* cf);

double GetValueFromGaussian(double x,double max,double zero);

void modify_pdb(char* infile, char* outfile, int exclude_het, int remove_water, int is_protein); // reorder protein atoms in PDB file
int rna_structure(char* infile);
int get_NextLine(char lines[][100], int nlines);
int is_rna_structure(char* infile);
int is_natural_amino(char* res);
int is_natural_nucleic(char* res);
void rewrite_residue2(char lines[][100], int nlines, int* wrote, FILE* outfile_ptr); // rewrite ordered residue

double get_yval(struct energy_matrix* energy_matrix, double area);

/*
  #ifdef __cplusplus
  }
  #endif
*/

#endif // include guard 
