#include "flexaid.h"
/***********************************************************************************
 * This subroutine returns the number of flexible dihedral bonds of an amino acid
 * given its three letter code in capital letters
 *
 ***********************************************************************************/

int number_of_dihedrals(char res[]){
  
  int ndihedrals;

  ndihedrals=0;

  if(strcmp(res,"ARG") == 0){ndihedrals=4;}
  if(strcmp(res,"ASN") == 0){ndihedrals=2;}
  if(strcmp(res,"ASP") == 0){ndihedrals=2;}
  if(strcmp(res,"CYS") == 0){ndihedrals=1;}
  if(strcmp(res,"GLN") == 0){ndihedrals=3;}
  if(strcmp(res,"GLU") == 0){ndihedrals=3;}
  if(strcmp(res,"HIS") == 0){ndihedrals=2;}
  if(strcmp(res,"ILE") == 0){ndihedrals=2;}
  if(strcmp(res,"LEU") == 0){ndihedrals=2;}
  if(strcmp(res,"LYS") == 0){ndihedrals=4;}
  if(strcmp(res,"MET") == 0){ndihedrals=3;}
  if(strcmp(res,"PHE") == 0){ndihedrals=2;}
  if(strcmp(res,"SER") == 0){ndihedrals=1;}
  if(strcmp(res,"THR") == 0){ndihedrals=1;}
  if(strcmp(res,"TRP") == 0){ndihedrals=2;}
  if(strcmp(res,"TYR") == 0){ndihedrals=2;}
  if(strcmp(res,"VAL") == 0){ndihedrals=1;}

  return ndihedrals;
}
