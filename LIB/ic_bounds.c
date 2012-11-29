#include "flexaid.h"

/****************************************************************************** 
 * SUBROUTINE ic_bounds calculates the bounds on the ic parameters for GPA to
 * restrict the search to a restricted area of the protein surface
 ******************************************************************************/
void ic_bounds(FA_Global* FA, char* rngopt){
	
	if(strcmp(rngopt,"loccen")==0 ||
	   strcmp(rngopt,"locclf")==0){
	  
		// set the number of spheres as bounds
		// so that every sphere can be selected
		// 0 is the INI position of the ligand
		FA->index_min = 1.0;
		FA->index_max = (double)FA->num_grd;

	}else if(strcmp(rngopt,"global")==0){
		FA->dis_min=FA->delta_angstron;
		FA->dis_max=FA->maxdst;
		FA->ang_min=-180.0;
		FA->ang_max=180.0;
		FA->dih_min=-180.0;
		FA->dih_max=180.0;      
	}
  

	if(FA->normal_modes != 0) {
		// 0 is the initial protein conformation
		FA->normalindex_min = 0.0;
		FA->normalindex_max = (double)FA->normal_grid_points;
	}


	return;
}
