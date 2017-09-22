#include "gaboom.h"
#include "maps.hpp"
#include "boinc.h"

/*****************************************************************/
/******* THIS FUNCTION IS USED TO SPLIT IN HALF THE   ************/
/******* SPACING BETWEEN 2 INTERSECTIONS OF THE GRID  ************/
/******* LEADING TO A MUCH MORE DENSE GRID THAT       ************/
/******* COVERS LESS VOLUME OF THE BINDING SITE       ************/
/*****************************************************************/

void slice_grid(FA_Global* FA,genlim* gene_lim,atom* atoms,resid* residue,gridpoint** cleftgrid) {

    std::map<std::string,int> cleftgrid_map;
    std::map<std::string,int>::iterator it,it2;
    std::map<std::string,int>::iterator begin,end;
	
	float sqrspa = pow(FA->spacer_length, 2.0f);
	float sqrhyp = pow(FA->spacer_length / sin(PI/4.0f), 2.0f);

	//printf("slicing grid...\n");
	
	for(int i=1; i<FA->num_grd; i++){
		std::string key = get_key((*cleftgrid)[i].coor);
        cleftgrid_map.insert(std::pair<std::string,int>(key, i));
	}
	
	begin = cleftgrid_map.begin();
	end = cleftgrid_map.end();
	
	for(it=begin; it!=end; ++it){
		for(it2=it; it2!=end; ++it2){
			
			if(it2 == it) continue;

			float sqrdis = sqrdist((*cleftgrid)[it->second].coor,(*cleftgrid)[it2->second].coor);
			if((fabs(sqrdis-sqrspa) < 0.001) ||
			   (fabs(sqrdis-sqrhyp) < 0.001)){
				
				float coor[3];
				coor[0] = ((*cleftgrid)[it->second].coor[0] + (*cleftgrid)[it2->second].coor[0]) / 2.0f;
				coor[1] = ((*cleftgrid)[it->second].coor[1] + (*cleftgrid)[it2->second].coor[1]) / 2.0f;
				coor[2] = ((*cleftgrid)[it->second].coor[2] + (*cleftgrid)[it2->second].coor[2]) / 2.0f;
				
				//string key1 = get_key((*cleftgrid)[it->second].coor);
				//string key2 = get_key((*cleftgrid)[it2->second].coor);
				std::string key = get_key(coor);

				/*
				printf("key1: %s\nkey2: %s\nkey: %s\n",
				       key1.c_str(), key2.c_str(), key.c_str());
				getchar();
				*/

				if(cleftgrid_map.find(key) == cleftgrid_map.end()){
					
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

                    cleftgrid_map.insert(std::pair<std::string,int>(key, FA->num_grd));

					FA->num_grd++;
				}
			}
		}
	}


	ic_bounds(FA,FA->rngopt);
	
        //Reset gene limit and gene length
        //Minimum is always 1.0 (not affected)
	FA->max_opt_par[0] = FA->index_max;
	
	gene_lim->max = FA->max_opt_par[0];
	set_bins(gene_lim);
	
	calc_cleftic(FA,(*cleftgrid));
	
        // Set new spacer_length
	FA->spacer_length /= 2.0f;
	
	return;
}