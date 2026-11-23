#include "flexaid.h"

/********************************************************************************
 * This function calculates the RSMD between atomic coordinates of the atoms in *
 * the register ori_ligatm and those for the atoms of the ligand in             *
 * residue[opt_res[0]] after reconstructing the coordinates using opt_par       *
 ********************************************************************************/

float calc_rmsp(int npar, float* icv, float* icu){
  float rmsp=0.0;
  int i;

  for(i=0;i<npar;i++){
    rmsp += (icv[i]-icu[i])*(icv[i]-icu[i]);
  }

  rmsp /= (float)npar;

  return  sqrt(rmsp);
}
