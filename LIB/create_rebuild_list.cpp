#include "flexaid.h"
#include "boinc.h"
/******************************************************************************
 * SUBROUTINE create_rebuild_list uses the map_par to create for each residue 
 * the list of atoms that will have to be reconstructed each time a ic2cf call 
 * is made.
 ******************************************************************************/
void create_rebuild_list(FA_Global* FA,atom* atoms,resid* residue){

	int i,j;
	int res=0;    /* residue to which an atom belongs */
	int newres;   /* temp variable for residue recontruction */
	int natm;

	for(i=0;i<FA->npar;i++){
        
        //printf("map_par[%d].typ=%d\n",i,FA->map_par[i].typ);
		newres=atoms[FA->map_par[i].atm].ofres;
        
        if(newres != res && residue[newres].type == 1){
            
			res=newres;
			FA->opt_res[FA->nors]=res;

			natm = residue[res].latm[0]-residue[res].fatm[0]+1;
            
			//printf("allocating memory for mov\n");
			FA->mov[FA->nors] = (int*)malloc(natm*sizeof(int));
			if(!FA->mov[FA->nors]){
				fprintf(stderr,"ERROR: memory allocation error for mov\n");
				Terminate(2);
			}
			memset(FA->mov[FA->nors],0,natm*sizeof(int));

			buildlist(FA,atoms,residue,res,0,&FA->nmov[FA->nors],FA->mov[FA->nors]);
            
            /*
			printf("nors=%d opt_res[%d]=%d nmov[%d]=%d\n",
			       FA->nors,FA->nors,FA->opt_res[FA->nors],FA->nors,FA->nmov[FA->nors]);
             */
        	if (FA->htpmode == false) {
        		for(j=0;j<FA->nmov[FA->nors];j++)
        		{
        			printf("mov[%d][%d]=%d\n",FA->nors,j,FA->mov[FA->nors][j]);
        		}
        	}
			//PAUSE;
			FA->nors++;
        }
	}
}