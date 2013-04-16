#include "flexaid.h"

/******************************************
 Assigns radius of atoms based on the
 atom types derived from SYBYL
 ******************************************/

void assign_radii_types(FA_Global* FA, atom* atoms, resid* residue)
{
    //  //const char Atoms_typename[11][5] = {"C3H0","C3H1","C4  ","N3H0","N3H1","N3H2","N4  ","O1H0","O2H1","S   ","DEF "};
    //const float radius[11] =              { 1.61f, 1.76f, 1.88f, 1.64f, 1.64f, 1.64f, 1.64f, 1.42f, 1.46f, 1.77f, 1.80f };

    for(int i=1; i<=FA->res_cnt; i++)
    {
        for(int j=residue[i].fatm[0]; j<=residue[i].latm[0]; j++)
        {
            switch(atoms[j].type){
                case 1:
                    break;
                case 2:
                    break;
                case 3:
                    break;
                case 4:
                    break;
                case 5:
                    break;
                case 6:
                    break;
                case 7:
                    break;
                case 8:
                    break;
                case 9:
                    break;
                case 10:
                    break;
                case 11:
                    break;
                case 12:
                    break;
                case 13:
                    break;
                case 14:
                    break;
                case 15:
                    break;
                case 16:
                    break;
                case 17:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.02f: 1.8f;
                    break;
                case 18:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.02f: 1.8f;
                    break;
                case 19:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.02f: 1.8f;
                    break;
                case 20:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.02f: 1.8f;
                    break;
                case 21:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.02f: 1.8f;
                    break;
                case 22:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.05f: 1.8f;
                    break;
                case 23:
                    atoms[j].radius = atoms[j].bonds[0] ? 0.64f: 1.47f;
                    break;
                case 24:
                    atoms[j].radius = atoms[j].bonds[0] ? 0.99f: 1.75f;
                    break;
                case 25:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.21f: 1.85f;
                    break;
                case 26:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.4f: 1.98f;
                    break;
                case 27:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.22f: 1.9f;
                    break;
                case 28:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.41f: 1.73f;
                    break;
                case 29:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.95f: 2.0f;
                    break;
                case 30:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.32f: 1.4f;
                    break;
                case 31:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.61f: 2.0f;
                    break;
                case 32:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.32f: 1.55f;
                    break;
                case 33:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.54f: 1.58f;
                    break;
                case 34:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.24f: 1.63f;
                    break;
                case 35:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.22f: 1.39f;
                    break;
                case 36:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.76f: 2.0f;
                    break;
                case 37:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.52f: 2.0f;
                    break;
                case 38:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.26f: 1.20f;
                    break;
                case 39:
                    atoms[j].radius = atoms[j].bonds[0] ? 1.5f: 2.0f;
                    break;
                case 40:
                    // solvent
                    break;
                
            }
                
        }
    }
}
