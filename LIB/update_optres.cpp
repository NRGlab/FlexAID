#include "flexaid.h"
#include "boinc.h"

/*****************************************************/
/******* this function updates the atoms str  ********/
/******* to assign the pointer pointing at    ********/
/******* the residue to optimize during dock  ********/
/*****************************************************/

void update_optres(atom* atoms, resid* residue, int atm_cnt, OptRes* optres_ptr,int num_optres)
{
    
	int i,j;
    
	for(i=0;i<num_optres;i++){
        
		for(j=1;j<=atm_cnt;j++){
			
			if(atoms[j].ofres == optres_ptr[i].rnum)
			{
				if(residue[atoms[j].ofres].type != 0 || !atoms[j].isbb)
				{
					atoms[j].optres = &optres_ptr[i];
					//printf("atom %d has optres\n", atoms[j].number);
				}
			}
            
		}
        
	}
    
	return;
    
}
