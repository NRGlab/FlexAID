#include "flexaid.h"

/******************************************
 Assigns radius of atoms based on the
 atom types derived from SYBYL
 ******************************************/
void assign_radii_types(FA_Global* FA, atom* atoms, resid* residue)
{
    for(int i=1; i<=FA->res_cnt; i++)
    {
        for(int j=residue[i].fatm[0]; j<=residue[i].latm[0]; j++)
        {
            switch(atoms[j].type){
                case 1: // C.1
                    atoms[j].radius = 1.88f;
                    break;
                case 2: // C.2
                    if((3 - atoms[j].bond[0]) > 0){
                        atoms[j].radius = 1.72f;
                    }else{
                        atoms[j].radius = 1.61f;
                    }
                    break;
                case 3: // C.3
                    atoms[j].radius = 1.88f;                    
                    break;
                case 4: // C.AR
                    atoms[j].radius = 1.76f;
                    break;
                case 5: // C.CAT
                    atoms[j].radius = 1.88f;                    
                    break;
                case 6: // N.1
                    atoms[j].radius = 1.64f;
                    break;
                case 7: // N.2
                    atoms[j].radius = 1.64f;
                    break;
                case 8: // N.3
                    atoms[j].radius = 1.64f;
                    break;
                case 9: // N.4
                    atoms[j].radius = 1.64f;
                    break;
                case 10: // N.AR
                    atoms[j].radius = 1.64f;
                    break;
                case 11: // N.AM
                    atoms[j].radius = 1.64f;
                    break;
                case 12: // N.PL3
                    atoms[j].radius = 1.64f;
                    break;
                case 13: // O.2
                    atoms[j].radius = 1.42f;
                    break;
                case 14: // O.3
                    atoms[j].radius = 1.46f;
                    break;
                case 15: // O.CO2
                    atoms[j].radius = 1.46f;
                    break;
                case 16: // O.AR
                    atoms[j].radius = 1.46f;
                    break;
                case 17: // S.2
                    atoms[j].radius = 1.782f;
                    break;
                case 18: // S.3
                    atoms[j].radius = 1.782f;
                    break;
                case 19: // S.O
                    atoms[j].radius = 1.782f;
                    break;
                case 20: // S.O2
                    atoms[j].radius = 1.782f;
                    break;
                case 21: // S.AR
                    atoms[j].radius = 1.782f;
                    break;
                case 22: // P.3
                    atoms[j].radius = 1.871f;
                    break;
                case 23: // F
                    atoms[j].radius = 1.560f;
                    break;
                case 24: // CL
                    atoms[j].radius = 1.735f;
                    break;
                case 25: // BR
                    atoms[j].radius = 1.978f;
                    break;
                case 26: // I
                    atoms[j].radius = 2.094f;
                    break;
                case 27: // SE
                    atoms[j].radius = 1.9f;                    
                    break;
                case 28: // MG
                    atoms[j].radius = 0.72f;
                    break;
                case 29: // SR
                    atoms[j].radius = 1.18f;
                    break;
                case 30: // CU
                    atoms[j].radius = 1.18f;
                    break;
                case 31: // MN
                    atoms[j].radius = 0.73f;
                    break;
                case 32: // HG
                    atoms[j].radius = 1.02f;
                    break;
                case 33: // CD
                    atoms[j].radius = 0.95f;
                    break;
                case 34: // NI
                    atoms[j].radius = 0.69f;
                    break;
                case 35: // ZN
                    atoms[j].radius = 0.74f;
                    break;
                case 36: // CA
                    atoms[j].radius = 1.00f;
                    break;
                case 37: // FE
                    atoms[j].radius = 0.61f;
                    break;
                case 38: // CO.OH
                    atoms[j].radius = 0.65f;
                    break;
                case 39: // DUMMY
                    atoms[j].radius = 2.0f;
                    break;
                case 40: // SOLVENT
                    break;
                
            }
            
            //printf("atom[%d].type=%d set to radius %.3f\n", atoms[j].number, atoms[j].type, atoms[j].radius);
        }
    }
}
