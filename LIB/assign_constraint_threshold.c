#include "flexaid.h"

//////////////////////////////////////////////////////////
///***       This function pre-calculates the       ***///
///***       maximal distances and angle threshold  ***///
///***       for constraints                        ***///
//////////////////////////////////////////////////////////


void assign_constraint_threshold(FA_Global* FA,atom* atoms,constraint* cons, int ncons)
{

  int i;

  float radA,radB;

  for(i=0;i<ncons;i++){

    cons[i].max_dist = 0.0f;
    cons[i].max_ang = 0.0f;

    radA = atoms[FA->num_atm[cons[i].anum1]].radius;
    radB = atoms[FA->num_atm[cons[i].anum2]].radius;

    if(radA != 0.0 && radB != 0.0)
      {

	cons[i].max_dist = 
	  (radA + radB +
	   2.0f*Rw);
	
	cons[i].max_ang =
	  (radA + radB +
	   2.0f*Rw);
	  
      }

  }

  return;

}
