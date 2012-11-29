#include "boinc.h"
#include "flexaid.h"

//////////////////////////////////////////////////////////
///// Function that builds the list of neighbouring //////
///// residues used to determine whether there is clash///
//////////////////////////////////////////////////////////

void build_close(FA_Global* FA, resid** residue, atom** atoms)
{

	int i,j;    // counters
	int kcalpha,calpha;   // C-alpha atm number
	int kres;     // source residue from which list is built

	int MIN_CLOSE_RESIDUES;
  
	// loop through list of flexible residues
	for(i=0;i<FA->nflxsc;i++){
    
		kres = FA->flex_res[i].inum;
		kcalpha = (*residue)[kres].fatm[0]+1;

		// all neighbours are rigid
		FA->flex_res[i].cflag = 0;


		MIN_CLOSE_RESIDUES = 10;

		FA->flex_res[i].close = (int*)malloc(MIN_CLOSE_RESIDUES*sizeof(int));
		if(!FA->flex_res[i].close){
			fprintf(stderr,"ERROR: memory allocation error for close\n");
			Terminate(2);
		}

		memset(FA->flex_res[i].close,0,5*sizeof(int));

		// loop through all protein residues
		for(j=1;j<=FA->res_cnt;j++){

			if(j == kres) continue;
    
			// skip HET groups
			if( (*residue)[j].type == 1 ) continue;

			// distance threshold is based on the Calpha atom only
			calpha = (*residue)[j].fatm[0]+1;

			if(dist((*atoms)[calpha].coor,(*atoms)[kcalpha].coor) <= MAX_CLOSE_DIST){
	
	
				if (FA->flex_res[i].nclose == MIN_CLOSE_RESIDUES) {

					//printf("re-allocating memory for close\n");
					MIN_CLOSE_RESIDUES += 5;
	  
					FA->flex_res[i].close = (int*)realloc(FA->flex_res[i].close,MIN_CLOSE_RESIDUES*sizeof(int));
					if(!FA->flex_res[i].close){
						fprintf(stderr,"ERROR: memory allocation error for close\n");
						Terminate(2);
					}
	
					memset(&FA->flex_res[i].close[MIN_CLOSE_RESIDUES-5],0,5*sizeof(int));

				}

				if ( (*atoms)[(*residue)[j].latm[0]].recs == 'm' ) { 
	  
					FA->flex_res[i].cflag = 1;

					FA->flex_res[i].close[FA->flex_res[i].nclose++] = j;
	  
				}else{

					// do not add rigid residues

				}
	
				/*
				  printf("found neighbour residue %s%c%d : latm.recs = %c\n",
				  (*residue)[j].name,(*residue)[j].chn,
				  (*residue)[j].number,(*atoms)[(*residue)[j].latm[0]].recs);
				*/

				continue;

			}
      

		}

    
		printf("Residue: %s %c %4d\n", (*residue)[kres].name,
		       (*residue)[kres].chn,(*residue)[kres].number);
		printf("CloseList (cflag=%d): ",FA->flex_res[i].cflag);
		for(j=0;j<FA->flex_res[i].nclose;++j){
			printf("\t%s %c %4d", 
			       (*residue)[FA->flex_res[i].close[j]].name,
			       (*residue)[FA->flex_res[i].close[j]].chn,
			       (*residue)[FA->flex_res[i].close[j]].number);    
		}
		printf("\n");
    
   
	}

	return;

}
