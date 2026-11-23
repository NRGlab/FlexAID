#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE assign_radius assigns the the atomic radius for all atoms present
 * in the system.
 ******************************************************************************/
float assign_radius(char atm[]){
  float r;
  char  c;

  c=atm[1];

  switch(c){

  case 'N':
    r=1.70f;
    break;
  case 'O':
    r=1.50f;
    break;
  case 'C':
    r=1.90f;
    break;
  case 'S':
    r=1.90f;
    break;
  case 'P':
    r=1.23f;
    break;
  case 'L':
    r=1.75f;
    break;
  case 'F':
    r=1.70f;
    break;
  case 'R':
    r=1.85f;
    break;
  case 'I':
    r=2.00f;
    break;
  case 'E':
    r=1.50f;
    break;
  case 'A':
    r=1.50f;
    break;
  default :
    r=1.50f;
    break;
  }

  return (r);
}
