#include "gaboom.h"
#include "maps.hpp"
#include "boinc.h"

/**************************************************************/
/******** THIS FUNCTION IS USED TO SELECT AREAS OF ************/
/******** THE GRID (INTERSECTIONS) THAT LEAD TO    ************/
/******** INDIVIDUALS WITH HIGHER CF VALUES        ************/
/******** AMONG THE TOP<POP_SELECTION> INDIVIDUALS ************/
/******** THIS FUNCTION ALSO EXPANDS THE GRID BY   ************/
/******** EXPANSION_FACTOR*SPACER THE CHOSEN INT.  ************/
/**************************************************************/

void partition_grid(FA_Global* FA,chromosome* chrom,genlim* gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid,int pop_selection,int expfac){

	std::map<std::string, std::vector<int> > cleftgrid_map;
	std::map<std::string, std::vector<int> >::iterator it;
	std::map<std::string, int > cleftgrid_struct_map;
	std::map<std::string, int >::iterator it2;
	std::vector<int>::iterator itv;

	//printf("partitionning grid\n");

	//Add grid intersection indexes for the <pop_selection> first chromosomes (based on evalue)
	//Store the indexes of selected chromosomes (key is the coordinate of the grid point)
	for(int i=0;i<pop_selection;i++){
		//Only consider the first gene
		//Find the integer value from 1 to num_grd
		int grd_idx = (int)chrom[i].genes->to_ic;

		std::string key = get_key((*cleftgrid)[grd_idx].coor);
		
		it = cleftgrid_map.find(key);
		if(it == cleftgrid_map.end()){
			std::vector<int> newvec;
			newvec.push_back(i);

			// key(coor) -> index individuals (ranked according to evalue)
			cleftgrid_map.insert(std::pair<std::string, std::vector<int> >(key, newvec));
		}else{
			it->second.push_back(i);
		}
	}

	//printf("partitionning grid: population inserted in map\n");

	int x,y,z;
	for(it=cleftgrid_map.begin(); it!=cleftgrid_map.end(); ++it){
		
		if(it->second.size() == 0) continue;

		for(x=-expfac;x<=expfac;x++){
			for(y=-expfac;y<=expfac;y++){
				for(z=-expfac;z<=expfac;z++){
					
					float coor[3];
					parse_key(it->first, coor);
					coor[0] += FA->spacer_length * (float)x;
					coor[1] += FA->spacer_length * (float)y;
					coor[2] += FA->spacer_length * (float)z;
					
					std::string key = get_key(coor);
					if(cleftgrid_map.find(key) == cleftgrid_map.end()){
						std::vector<int> emptyvec;
					        cleftgrid_map.insert(std::pair<std::string, std::vector<int> >(key, emptyvec));
					}
				}
			}
		}
	}
	
	//printf("partitionning grid: expanded partition\n");

	FA->num_grd = 1;

	// all keys represent a unique grid point in the new partitionned grid map
	// reinsert each key into cleftgrid structure (start at index = 1, 0 is reference)
	// increase size of cleftgrid structure if necessary
	for(it=cleftgrid_map.begin(); it!=cleftgrid_map.end(); ++it){
		float coor[3];
		parse_key(it->first,coor);
		
		if (FA->num_grd==FA->MIN_CLEFTGRID_POINTS){
			FA->MIN_CLEFTGRID_POINTS *= 2;
								
			(*cleftgrid) = (gridpoint*)realloc((*cleftgrid),FA->MIN_CLEFTGRID_POINTS*sizeof(gridpoint));
			if ((*cleftgrid) == NULL){
				fprintf(stderr,"ERROR: memory reallocation error for cleftgrid (partition)\n");
				Terminate(2);
			}
		}		

		(*cleftgrid)[FA->num_grd].coor[0] = coor[0];
		(*cleftgrid)[FA->num_grd].coor[1] = coor[1];
		(*cleftgrid)[FA->num_grd].coor[2] = coor[2];
		
		// key(coor) -> grid index in structure
		cleftgrid_struct_map.insert(std::pair<std::string,int>(it->first, FA->num_grd));
		
		FA->num_grd++;
	}
  
	//printf("partitionning grid: insert into structure done\n");

	ic_bounds(FA,FA->rngopt);
  
	//Reset gene limit and gene length
	//Minimum is always 1.0 (not affected)
	FA->max_opt_par[0] = FA->index_max;

	gene_lim->max = FA->max_opt_par[0];
	set_bins(gene_lim);

	// adjust population of chromosomes
	// loop through grid points in original grid
	for(it=cleftgrid_map.begin(); it!=cleftgrid_map.end(); ++it){
		// loop through individuals of the given grid point
		for(itv=it->second.begin(); itv!=it->second.end(); itv++){
			// what is the index of the given grid point in the new cleftgrid structure
			if((it2 = cleftgrid_struct_map.find(it->first)) != cleftgrid_struct_map.end()){
				
				chrom[*itv].genes->to_ic = (double)it2->second;
				chrom[*itv].genes->to_int32 = ictogene(gene_lim,it2->second);

			}else{
				fprintf(stderr, "ERROR: Could not find key in cleftgrid_struct_map\n");
				Terminate(22);
			}
		}
	}
  
	calc_cleftic(FA,(*cleftgrid));

}