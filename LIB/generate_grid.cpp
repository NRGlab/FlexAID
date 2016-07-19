#include "flexaid.h"
#include "boinc.h"
#include "maps.hpp"

/*******************************************/
/*  generates a grid from the spheres of 
    loccen or spheres combination from
    locclf (using a cleft)                 */
/*******************************************/

gridpoint* generate_grid(FA_Global* FA,sphere* spheres, atom* atoms, resid* residue){
	
	float sqrrad;
	float c[3], min[3], max[3];
	gridpoint* cleftgrid = NULL;
	std::map<std::string,int> cleftgrid_map;

	cleftgrid = (gridpoint*)malloc(FA->MIN_CLEFTGRID_POINTS*sizeof(gridpoint));    
	if (cleftgrid == NULL){
		fprintf(stderr,"ERROR: memory allocation error for cleftgrid\n");
		Terminate(2);
	}

    cleftgrid[0].number = 0;
    for (int j=0;j<3;j++) cleftgrid[0].coor[j] = atoms[residue[FA->res_cnt].gpa[0]].coor[j];
    cleftgrid[0].dis = atoms[residue[FA->res_cnt].gpa[0]].dis;
    cleftgrid[0].ang = atoms[residue[FA->res_cnt].gpa[0]].ang;
    cleftgrid[0].dih = atoms[residue[FA->res_cnt].gpa[0]].dih;
    
	printf("will build a grid with spacing %.3f\n", FA->spacer_length);

	FA->num_grd = 1; // set counter to 1 because 0 is the INI conformation of the ligand
	while(spheres != NULL){
		
        /*
		cout << "CENTER=" << spheres->center[0] << " " << spheres->center[1] << " " << spheres->center[2] << endl;
		cout << "RADIUS=" << spheres->radius << endl;
		*/
        
		if ( (float)( 1.0 / FA->spacer_length ) - (float)((int)( 1.0 / FA->spacer_length )) > 0.001 ){
			min[0] = (float)((int)( (spheres->center[0] - spheres->radius) / FA->spacer_length )) * FA->spacer_length;
			min[1] = (float)((int)( (spheres->center[1] - spheres->radius) / FA->spacer_length )) * FA->spacer_length;
			min[2] = (float)((int)( (spheres->center[2] - spheres->radius) / FA->spacer_length )) * FA->spacer_length;
			max[0] = (float)((int)( (spheres->center[0] + spheres->radius) / FA->spacer_length ) + 1.0) * FA->spacer_length;
			max[1] = (float)((int)( (spheres->center[1] + spheres->radius) / FA->spacer_length ) + 1.0) * FA->spacer_length;
			max[2] = (float)((int)( (spheres->center[2] + spheres->radius) / FA->spacer_length ) + 1.0) * FA->spacer_length;
		}else{
			min[0] = (float)((int)( (spheres->center[0] - spheres->radius - FA->spacer_length )));
			min[1] = (float)((int)( (spheres->center[1] - spheres->radius - FA->spacer_length )));
			min[2] = (float)((int)( (spheres->center[2] - spheres->radius - FA->spacer_length )));
			max[0] = (float)((int)( (spheres->center[0] + spheres->radius + FA->spacer_length ) + 1.0));
			max[1] = (float)((int)( (spheres->center[1] + spheres->radius + FA->spacer_length ) + 1.0));
			max[2] = (float)((int)( (spheres->center[2] + spheres->radius + FA->spacer_length ) + 1.0));
		}

		c[0] = min[0];
		c[1] = min[1];
		c[2] = min[2];
		
		sqrrad = spheres->radius * spheres->radius;

		while(c[2] < max[2]){
			while(c[1] < max[1]){
				while(c[0] < max[0]){

					if(sqrdist(spheres->center,c) < sqrrad){
						
						std::string key = get_key(c);
						
						if(cleftgrid_map.find(key) == cleftgrid_map.end()){							
							if (FA->num_grd==FA->MIN_CLEFTGRID_POINTS){
								FA->MIN_CLEFTGRID_POINTS *= 2;
								
								cleftgrid = (gridpoint*)realloc(cleftgrid,FA->MIN_CLEFTGRID_POINTS*sizeof(gridpoint));
								if (cleftgrid == NULL){
									fprintf(stderr,"ERROR: memory reallocation error for cleftgrid\n");
									Terminate(2);
								}
							}		
							
							cleftgrid[FA->num_grd].coor[0] = c[0];
							cleftgrid[FA->num_grd].coor[1] = c[1];
							cleftgrid[FA->num_grd].coor[2] = c[2];
							
							FA->num_grd++;
                            cleftgrid_map.insert(std::pair<std::string,int>(key, FA->num_grd));

						}

					}
					
					c[0] += FA->spacer_length;
				}

				c[0] = min[0];
				c[1] += FA->spacer_length;
			}
			
			c[0] = min[0];
			c[1] = min[1];
			c[2] += FA->spacer_length;
			
		}
		
		spheres = spheres->prev;
	}

	printf("built a grid with %d vertices\n", FA->num_grd - 1);

	return cleftgrid;
}
