#include "flexaid.h"

/*********************************************/
/* This function alters the backbone of the  */
/* protein. It uses a row in the grid to     */
/* select modes and amplitudes of vectors    */
/* The function takes for arguments the      */
/* initial protein conformation              */
/*********************************************/

void alter_mode(atom* atoms,resid* residue, float *grid, int res_cnt, int nmodes) {
	int i,j,k,l;
	int rot;

	float amplitude;
  
	for (i=1;i<=res_cnt;++i) {
		if(residue[i].type == 1){ continue; }
    
		rot=residue[i].rot;
		for(j=residue[i].fatm[rot];j<=residue[i].latm[rot];++j){

			// reset original coordinates
			for(l=0;l<3;l++){ atoms[j].coor[l] = atoms[j].coor_ori[l]; }

			for (k=0;k<nmodes;++k) {
				amplitude = grid[k];
				//printf("For atom[%d] Ampli:%f\n",j,amplitude);
	
				for(l=0;l<3;l++){ atoms[j].coor[l] += amplitude * atoms[j].eigen[k][l]; }
			} 
		}
	}
  
	return;
}

