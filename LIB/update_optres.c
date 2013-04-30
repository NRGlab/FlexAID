#include "flexaid.h"
#include "boinc.h"

/*****************************************************/
/******* this function updates the atoms str  ********/
/******* to assign the pointer pointing at    ********/
/******* the residue to optimize during dock  ********/
/*****************************************************/

void update_optres(atom* atoms, int atm_cnt, OptRes* optres_ptr,int num_optres)
{
    
    int i,j;
    
    for(i=0;i<num_optres;i++){
        
        for(j=1;j<=atm_cnt;j++){
            
            if(atoms[j].ofres == optres_ptr[i].rnum)
            {
                atoms[j].optres = &optres_ptr[i];
            }
            
        }
        
    }
    
    return;
    
}
