#include "flexaid.h"

/***************************************************************************************
 * This function checks if there are steric clashes between the atoms in the list and
 * other atoms in the protein which are rigid (i.e., not rotamers and ligand)
 **************************************************************************************/
int check_clash(FA_Global* FA,atom* atoms,resid* residue,int res_cnt,int total, int list[]){
	int i,j,k;
	float dis,rad;
	//  float frac=1.0;

	for(k=0;k<total;k++){

		/*
		printf("==============================\n");
		printf("FOR ATOM %s %c %d %s\n", 
		       residue[atoms[list[k]].ofres].name,
		       residue[atoms[list[k]].ofres].chn,
		       residue[atoms[list[k]].ofres].number,
		       atoms[list[k]].name);
		*/

		for(i=1;i<=res_cnt;i++){
    
			for(j=residue[i].fatm[0];j<=residue[i].latm[0];j++){
      
				if(atoms[j].recs == 'r'){
					
					// skip internal clash
					if(atoms[list[k]].ofres != i){

						dis=distance(atoms[j].coor,atoms[list[k]].coor);
						rad=atoms[j].radius+atoms[list[k]].radius;
						//if(dis<=2.0){
	    
						
						if(dis < FA->rotamer_permeability*rad){
	      
							/*
							printf("RES %4d %s\n",residue[i].number,residue[i].name);
							printf("RES %4d %s ATM %5d %4s %c %f %f %f\n",
							       residue[atoms[j].ofres].number,residue[atoms[j].ofres].name,
							       atoms[j].number,atoms[j].name,atoms[j].recs,
							       atoms[j].coor[0],atoms[j].coor[1],atoms[j].coor[2]);
	      
							printf("RES %4d %s ATM %5d %4s %c %f %f %f\n",
							       residue[atoms[list[k]].ofres].number,residue[atoms[list[k]].ofres].name,
							       atoms[list[k]].number,atoms[list[k]].name,atoms[list[k]].recs,
							       atoms[list[k]].coor[0],atoms[list[k]].coor[1],atoms[list[k]].coor[2]);

							printf("DIS %8.4f RAD %8.4f DEE_CLA %8.4f\n",dis,rad,rad*FA->dee_clash);
							*/
	      	      
							if ( FA->useflexdee && dis > FA->dee_clash*rad ) {
								return (0);
							} else {
								return (1);
							}

						}

					}

				}

			}
		}
	}

	return (0);
}
