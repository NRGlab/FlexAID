#include <sstream>
#include <string>
#include <vector>

#include "flexaid.h"
#include "boinc.h"

/****************************************/
/* assign_flexclash: assigns the flexible
   bonds that are within the shortest path
   between 2 atoms. this function is
   necessary for the DEE in the ligand
   flexibility                          */
/****************************************/


std::vector<std::string> split_string(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream token_stream(str);

    while (std::getline(token_stream, token, delimiter)) {
        tokens.push_back(token);
    }

    return tokens;
}

void assign_shortflex(resid* residue, int tot, int fdih, atom* atoms)
{
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
			std::vector<std::string> iatmnum;
			std::vector<std::string>::iterator it,it2;

			iatmnum = split_string(residue->shortpath[i][j], ' ');

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
}