/* calculates the internal coordinates of the atom sent to parameter list. References points are used */

#include "flexaid.h"

void calc_cleftic(FA_Global* FA,gridpoint* cleftgrid){

  float ci[3],ca[3],cb[3],cc[3];
  int i,j; // counters

  ca[0] = FA->ori[0]+1;
  ca[1] = FA->ori[1];
  ca[2] = FA->ori[2];

  cb[0] = FA->ori[0];
  cb[1] = FA->ori[1];
  cb[2] = FA->ori[2];
  
  cc[0] = FA->ori[0];
  cc[1] = FA->ori[1]+1;
  cc[2] = FA->ori[2];

  for (i=1;i<FA->num_grd;i++){

    for (j=0;j<3;j++){
      ci[j] = cleftgrid[i].coor[j];
    }
    
    cleftgrid[i].dis = distance(ci,ca);
    cleftgrid[i].ang = bndang(ci,ca,cb);
    cleftgrid[i].dih = dihedral(ci,ca,cb,cc);
  }
  
  return;
}
