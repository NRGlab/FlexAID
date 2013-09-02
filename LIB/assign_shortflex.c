#include "boost/algorithm/string.hpp"

#include "flexaid.h"
#include "boinc.h"

using namespace std;

/****************************************/
/* assign_flexclash: assigns the flexible
   bonds that are within the shortest path
   between 2 atoms. this function is
   necessary for the DEE in the ligand
   flexibility                          */
/****************************************/

void assign_shortflex(resid* residue, int tot, int fdih, atom* atoms)
{
	int fatm = residue->fatm[0];
	
	residue->shortflex = (int***)malloc(tot*sizeof(int**));
	if(!residue->shortflex){
		fprintf(stderr,"ERROR: Could not allocate memory for residue->shortflex\n");
		Terminate(2);
	}else{
		for(int i=0; i<tot; i++){
			residue->shortflex[i] = (int**)malloc(tot*sizeof(int*));
			if(!residue->shortflex[i]){
				fprintf(stderr,"ERROR: Could not allocate memory for residue->shortflex[%d]\n", i);
				Terminate(2);
			}else{
				for(int j=0; j<tot; j++){
					residue->shortflex[i][j] = (int*)malloc((fdih+1)*sizeof(int));
					if(!residue->shortflex[i][j]){
						fprintf(stderr,
							"ERROR: Could not allocate memory for residue->shortflex[%d][%d]\n", i, j);
						Terminate(2);
					}else{
						memset(residue->shortflex[i][j],-1,(fdih+1)*sizeof(int));
					}
				}
			}
		}
	}
	
	for(int i=0; i<tot; i++){
		for(int j=0; j<tot; j++){
			vector<string> iatmnum;
			vector<string>::iterator it,it2;
			boost::split(iatmnum, residue->shortpath[i][j], boost::is_any_of(" "));
			
			//cout << residue->shortpath[i][j] << endl;
			
			// atoms must at least be seperated by a dihedral
			if(iatmnum.size() > 3){
				int counter = 0;
				for(it=iatmnum.begin(); it!=iatmnum.end(); ++it){
					it2 = it; ++it2;
					if(it!=iatmnum.begin() && it2!=iatmnum.end() && *it2!=iatmnum.back()){
						int atm1 = atoi((*it).c_str());
						int atm2 = atoi((*it2).c_str());
						
						// loops through flexible bonds
						for(int k=1; k<=fdih; k++){
							int flexatm = residue->bond[k];
							
							if((atoms[flexatm].rec[0] == atm1 && atoms[flexatm].rec[1] == atm2) ||
							   (atoms[flexatm].rec[0] == atm2 && atoms[flexatm].rec[1] == atm1)){
							/*
							if((atoms[atoms[flexatm].rec[0]].number == atm1 && 
							    atoms[atoms[flexatm].rec[1]].number == atm2) ||
							   (atoms[atoms[flexatm].rec[0]].number == atm2 && 
							    atoms[atoms[flexatm].rec[1]].number == atm1)){
							*/
								/*
								cout << "added " << k << " "
								     << atoms[atm1].number << " "
								     << atoms[atm2].number << endl;
								*/
								
								residue->shortflex[i][j][counter++] = k;
								break;
							}
						}
					}
				}
			}
		}
	}
	
	return;
}
