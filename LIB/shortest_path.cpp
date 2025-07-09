#include "flexaid.h"
#include "boinc.h"

/*******************************************/
/* Shortest_path: finds the shortest path  */
/* from any atom in residue to all other   */
/* atoms of the same residue               */
/*******************************************/

void shortest_path(resid* residue, int tot, atom* atoms)
{	
	int fatm = residue->fatm[0];
	//int latm = residue->latm[0];
	
	if(residue->shortpath == NULL){
		residue->shortpath = (char***)malloc(tot*sizeof(char**));
		if(!residue->shortpath){
			fprintf(stderr,"ERROR: Could not allocate memory for residue->shortpath\n");
			Terminate(2);
		}
	}
	
	for(int i=0; i<tot; i++){
		residue->shortpath[i] = (char**)malloc(tot*sizeof(char*));
		if(!residue->shortpath[i]){
			fprintf(stderr,"ERROR: Could not allocate memory for residue->shortpath[%d]\n", i);
			Terminate(2);
		}else{
			for(int j=0; j<tot; j++){
				// 90001-90002
				residue->shortpath[i][j] = (char*)malloc(MAX_SHORTEST_PATH*6*sizeof(char));
				if(!residue->shortpath[i][j]){
					fprintf(stderr,
						"ERROR: Could not allocate memory for residue->shortpath[%d][%d]\n", i,j);
					Terminate(2);
				}
			}
		}
	}

	for(int i=0; i<tot; i++){
		//cout << "objective " << atoms[fatm+i].number << endl;
		for(int j=0; j<tot; j++){
			
			std::stringstream ss("");
			
			if(j==i){
				//ss << atoms[fatm+i].number;
				ss << fatm+i;
				strcpy(residue->shortpath[i][j], ss.str().c_str());
				/*
				  printf("shortest path from %d to %d is %s\n",
				  atoms[fatm+j].number, atoms[fatm+i].number, residue->shortpath[i][j]);
				*/
				continue;
			}
			
			std::vector<int> queue;
			std::map<int, std::string> mark;
			std::map<int, std::string>::iterator it;
			queue.push_back(j);
			
			//ss << atoms[fatm+j].number;
			ss << fatm+j;
			mark[j] = ss.str();
			
			while(queue.size()){
				int qa = queue.front();
				queue.erase(queue.begin());
				
				//cout << "dequeued " << atoms[fatm+qa].number << endl;
				if(qa == i){
					//ss.str("");
					//ss << mark[qa] << ' ' << fatm+i;
					//ss << mark[qa] << ' ' << atoms[fatm+i].number;
					//mark[qa] = ss.str();
					
					if(mark[qa].length() > MAX_SHORTEST_PATH*6){
						fprintf(stderr,"ERROR: Shortest path buffer too short for the size of the molecule\n");
						Terminate(25);
					}
					strcpy(residue->shortpath[j][qa],mark[qa].c_str());
					/*
					  printf("shortest path from %d to %d is %s\n",
					  atoms[fatm+j].number, atoms[fatm+qa].number, mark[qa].c_str());
					*/
					break;
				}
				
				for(int k=0; k<tot; k++){
					// add all atoms bonded to atom j 
					/*
					  printf("bonded[%d][%d]=%d\n", 
					  atoms[fatm+qa].number, 
					  atoms[fatm+k].number,
					  residue->bonded[qa][k]);
					*/
					if(residue->bonded[qa][k] == 1){
						// unless atom has already been visited (marked)
						it = mark.find(k);
						if(it == mark.end()){
							queue.push_back(k);
							//cout << "enqueued " << atoms[fatm+k].number << endl;
							
							ss.str("");
							//ss << mark[qa] << ' ' << atoms[fatm+k].number;
							ss << mark[qa] << ' ' << fatm+k;
							mark[k] = ss.str();
						}
					}
				}
			}
		}
	}	
}

/* BFS algorithm
   1  procedure BFS(G,v):
   2      create a queue Q
   3      enqueue v onto Q
   4      mark v
   5      while Q is not empty:
   6          t ← Q.dequeue()
   7          if t is what we are looking for:
   8              return t
   9          for all edges e in G.adjacentEdges(t) do
   10             u ← G.adjacentVertex(t,e)
   11             if u is not marked:
   12                  mark u
   13                  enqueue u onto Q
   14     return none
*/
