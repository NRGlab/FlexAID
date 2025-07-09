#include "flexaid.h"
#include "boinc.h"
/***********************************************************************************************
 * This function build for each residue its rotamers.
 *
 **********************************************************************************************/
void build_rotamers(FA_Global* FA,atom** atoms,resid* residue,rot* rotamer){
	int   i,j,k,l,m,j1,n;
	int   kres;
	int   buildcc_list[15];

	int   total_build;
	int   delta_atm;
	int   flag;
	int   rigid_clash=0;
	int   old_atm;

	int   old_num_ref=0;
	int   natm=0;

	int   MIN_ROTAMER=0;
	int   MIN_NEWROT=10;

	int   nbonded=0;
	int   bondlist[MAX_ATM_HET];
	int   neighbours[MAX_ATM_HET];
  
	int   nindex;

	int   num_ref=0;       // builds a PDBnum to map rotamer atoms in num_atm

	FILE* outrot = NULL;
					

	// Create a PDB number used for new rotamer atoms
	//num_ref=(*atoms)[residue[FA->res_cnt].latm[residue[FA->res_cnt].trot]].number;
  
	FA->nflxsc_real = 0 ;
  
	// set side-chain atoms as mobile for flexible residues
	for(k=0;k<FA->nflxsc;k++){
		kres=FA->flex_res[k].inum;
    
		// find C-beta atom (CB is always rigid)
		flag=0;
		for(i=residue[kres].fatm[0];i<=residue[kres].latm[0];i++){
			if(flag){
				(*atoms)[i].recs='m';
				
				/*
				printf("atom '%s' from %s %c %4d set as mobile (dih=%.3f)\n",
				       (*atoms)[i].name,
				       residue[kres].name,residue[kres].chn,residue[kres].number,
				       (*atoms)[i].dih);
				*/
			}
			
			if(strcmp(" CB ",(*atoms)[i].name) == 0){flag=1;}
		}
	}  
	
	//printf("Testing rotamers for clash... (libsize=%d)\n",FA->rotlibsize);
	
	for(i=0;i<FA->nflxsc;i++){
		//printf("---------\nFLEXSC(%d): %d %c %s\n",i,FA->flex_res[i].num,FA->flex_res[i].chn,FA->flex_res[i].name);
		MIN_ROTAMER = FA->MIN_ROTAMER;

		kres=FA->flex_res[i].inum;
		
		buildic(FA,*atoms,residue,kres);
		
		/*
		  printf("RESIDUE(%d): %d %c %s\n",k,residue[kres].number,residue[kres].chn,residue[kres].name);
		  printf("fatm=%d latm=%d fdih=%d\n",residue[kres].fatm[0],residue[kres].latm[0],residue[kres].fdih);
		  printf("TOTAL NUMBER OF ATOMS: %d\n",FA->atm_cnt);
		*/

		//for(j=1;j<=residue[kres].fdih;j++){
		//printf("Bond(%d)=%d %f ",j,residue[kres].bond[j],(*atoms)[residue[kres].bond[j]].dih);
		//}
		//printf("\n");
		//PAUSE;
		l=0;
		n=0;
		nindex=0;
		while(l<FA->rotlibsize){
			
			if(strcmp(rotamer[l].res,residue[kres].name) == 0){
				
				atom* at_cb = NULL;

				/*
				printf("---------------------\n");
				printf("Rotamer %d for %s %d %c \n",rotamer[l].nid,rotamer[l].res,
				       residue[kres].number,residue[kres].chn);	
				*/
				
				residue[kres].trot++;
	
				if(residue[kres].trot == MIN_ROTAMER){
					//printf("allocating memory for fatm/latm (rotamer)\n");
			  
					MIN_ROTAMER+=MIN_NEWROT;

					residue[kres].fatm = (int*)realloc(residue[kres].fatm,MIN_ROTAMER*sizeof(int));
					residue[kres].latm = (int*)realloc(residue[kres].latm,MIN_ROTAMER*sizeof(int));
		
					if(!residue[kres].fatm || !residue[kres].latm){
						fprintf(stderr,"ERROR: Could not re-allocate memory for residue[%d].fatm/latm\n",kres);
						Terminate(2);
					}
	  
					memset(&residue[kres].fatm[MIN_ROTAMER-MIN_NEWROT],0,MIN_NEWROT*sizeof(int));
					memset(&residue[kres].latm[MIN_ROTAMER-MIN_NEWROT],0,MIN_NEWROT*sizeof(int));
				}
	
				residue[kres].rot = residue[kres].trot;
				
				total_build=0;
				delta_atm=FA->atm_cnt-residue[kres].fatm[0]+1;
				
				residue[kres].fatm[residue[kres].trot]=residue[kres].fatm[0]+delta_atm;
				residue[kres].latm[residue[kres].trot]=residue[kres].latm[0]+delta_atm;
				
				old_num_ref=num_ref;
				old_atm=FA->atm_cnt;
				
				for(j=residue[kres].fatm[0];j<=residue[kres].latm[0];j++){
					FA->atm_cnt++;
	  
					//printf("new atom from %d (atm_cnt=%d)\n",j,FA->atm_cnt);
	  
					if(FA->atm_cnt == FA->MIN_NUM_ATOM){
	    
						FA->MIN_NUM_ATOM += 500;
						//printf("new MIN_NUM_ATOM: %d\n", FA->MIN_NUM_ATOM);

						(*atoms) = (atom*)realloc((*atoms),FA->MIN_NUM_ATOM*sizeof(atom));

						if(!(*atoms)){
							fprintf(stderr,"ERROR: Could not allocate memory for atoms\n");
							Terminate(2);
						}
	    
						memset(&(*atoms)[FA->MIN_NUM_ATOM-500],0,500*sizeof(atom));
					}
	  
	  
					num_ref++;
	  
					if((*atoms)[j].recs == 'm'){
						buildcc_list[total_build++]=j+delta_atm;
	    
						if((*atoms)[j].rec[3] != 0){
							(*atoms)[j].shift = (*atoms)[(*atoms)[j].rec[3]].dih - (*atoms)[j].dih;
						}
					}else if(!strcmp((*atoms)[j].name," CB ")){ at_cb = &(*atoms)[j]; }
	  
					(*atoms)[j+delta_atm] = (*atoms)[j];
					//(*atoms)[j+delta_atm].number += delta_atm;
					//(*atoms)[j+delta_atm].number = num_ref;
	  
					//FA->num_atm[num_ref] = j+delta_atm;
					
					for(m=0;m<4;m++){
						if((*atoms)[j+delta_atm].rec[m] != 0){
							(*atoms)[j+delta_atm].rec[m] += delta_atm;
						}
					}
				}
	
	
				for(j=1;j<=residue[kres].fdih;j++){
	  
					j1=residue[kres].bond[j]+delta_atm;
	  
					(*atoms)[j1].dih=rotamer[l].chi[j-1];
	  
					if((*atoms)[j1].dih > 180.0){
						(*atoms)[j1].dih -= 360.0;
					}
					if((*atoms)[j1].dih < -180.0){
						(*atoms)[j1].dih += 360.0;
					}
	  
					(*atoms)[(*atoms)[j1].rec[3]].dih=(*atoms)[j1].dih+(*atoms)[j1].shift;
	  
					if((*atoms)[(*atoms)[j1].rec[3]].dih > 180.0){
						(*atoms)[(*atoms)[j1].rec[3]].dih -= 360.0;
					}
	  
					if((*atoms)[(*atoms)[j1].rec[3]].dih < -180.0){
						(*atoms)[(*atoms)[j1].rec[3]].dih += 360.0;
					}
	  
				}

				/*
				printf("\nANGLES: ");
				for(j=1;j<=residue[kres].fdih;j++){
					j1=residue[kres].bond[j]+delta_atm;
					printf("%8.3f ",(*atoms)[j1].dih);
				}
				printf("\n");
				*/

				buildcc(FA,*atoms,total_build,buildcc_list);
				
				rigid_clash=check_clash(FA,(*atoms),residue,FA->res_cnt,total_build,buildcc_list);
				
				if(rigid_clash == 1){
					
					//printf("REJECTED DUE TO STERIC CLASHES WITH RIGID ATOM(S)\n\n");
					//getchar();	
					
					residue[kres].trot--;
					FA->atm_cnt=old_atm;
					num_ref=old_num_ref;
	  
				}else{
					
					//printf("ACCEPTED\n\n");
					
					if(FA->rotout){
						if(outrot == NULL){
							outrot = fopen("rotamers.pdb", "w");
							if(!outrot){
								fprintf(stderr, "ERROR: Could not open 'rotamers.pdb' for writing\n");
								outrot = stderr;
							}
						}
						
						fprintf(outrot, "MODEL       %2d\n", residue[kres].trot);
						resid* res = &residue[at_cb->ofres];
						fprintf(outrot, "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n",
							at_cb->number, at_cb->name,
							res->name, res->chn, res->number,
							at_cb->coor[0], at_cb->coor[1], at_cb->coor[2] );
						for(int m=0; m<total_build; m++){
							atom* at = &(*atoms)[buildcc_list[m]];
							fprintf(outrot, "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f\n",
								at->number, at->name,
								res->name, res->chn, res->number,
								at->coor[0], at->coor[1], at->coor[2] );
						}
						fprintf(outrot, "ENDMDL\n");
					}
					
					if(FA->nrg_suite){
						printf("Rotamer for %s%d%c with dihedrals", 
						       residue[kres].name, residue[kres].number,
						       residue[kres].chn == ' ' ? '-' : residue[kres].chn);
						
						for(j=0;j<residue[kres].fdih;j++){
							printf(" %8.3f", rotamer[l].chi[j]);
						}
						printf("\n");
					}
				}
				
				residue[kres].rot = 0;
				
				n++;
				
			}
      
			l++;
      
		}
    
		printf("%d possible rotamer(s) for residue %s %d %c\n",
			   residue[FA->flex_res[i].inum].trot, residue[FA->flex_res[i].inum].name,
			   residue[FA->flex_res[i].inum].number, 
			   residue[FA->flex_res[i].inum].chn == ' ' ? '-' : residue[FA->flex_res[i].inum].chn);
		
		// if no rotamers were accepted , set as rigid all atoms of side-chain
		if (residue[FA->flex_res[i].inum].trot == 0) {
			
			for(j=residue[FA->flex_res[i].inum].fatm[0];j<=residue[FA->flex_res[i].inum].latm[0];j++){
				
				(*atoms)[j].recs = 'r';
				
			}
			
		}else{
			
			FA->nflxsc_real++;
			
			/*
			  for(k=0; k<=residue[kres].trot; k++) { 
			  printf("rot[%d] fatm=%d latm=%d\n",
			  k,residue[kres].fatm[k],residue[kres].latm[k]);
			  }
			*/

		}

		fflush(stdout);
		
	}
	
	
	// set side-chain atoms as rigid for side-chains in which no rotamers were accepted
	// for those in which at least one rotamer was accepted
	// add to optimized residue struct and build bondedlist accordingly
	for(k=0;k<FA->nflxsc;k++){
		kres=FA->flex_res[k].inum;

		if(residue[kres].trot > 0 && residue[kres].bonded == NULL){
      
			if(!FA->score_ligand_only){
				FA->optres = (OptRes*)realloc(FA->optres,FA->MIN_OPTRES*sizeof(OptRes));
				if(!FA->optres){
					fprintf(stderr,"ERROR: Could not re-allocate memory for FA->optres.\n");
					Terminate(2);
				}
				memset(&FA->optres[FA->MIN_OPTRES-1],0,sizeof(OptRes));
				
				natm = residue[kres].latm[0]-residue[kres].fatm[0]+1;
				
				
				//printf("residue[%d].fatm = %d\n",residue[kres].number,residue[kres].fatm[0]);
				
				FA->optres[FA->MIN_OPTRES-1].rnum = kres;
				// 0: side-chain (protein)
				FA->optres[FA->MIN_OPTRES-1].type = 0;
				FA->optres[FA->MIN_OPTRES-1].tot = natm;
				
				for(i=0;i<natm;i++){
					
					nbonded=0;
					
					bondedlist(*atoms,residue[kres].fatm[0]+i,FA->bloops,&nbonded,bondlist,neighbours);
					
					/*
					  printf("bondlist[%d] for atom %d\t",i,(*atoms)[residue[kres].fatm[0]+i].number);
					  for(j=0;j<nbonded;j++) printf("%6d(%2d)",(*atoms)[bondlist[j]].number,neighbours[j]);
					  printf("\n");
					*/
					
					update_bonded(&residue[kres],natm,nbonded,bondlist,neighbours);
				}
				
				// prints bonded matrix
				/*
				  for(i=0;i<natm;i++){
				  printf("\t%5d\t",i);
				  for(j=0;j<natm;j++){
				  printf("%2d",residue[kres].bonded[i][j]);
				  }
				  printf("\n");
				  }
				  getchar();
				*/
				
				FA->num_optres++;
				
				FA->MIN_OPTRES++;
			}
		}else{
    
			for(i=residue[kres].fatm[0];i<=residue[kres].latm[0];i++){
				(*atoms)[i].recs='r';
			}
      
		}
    
	}  

	if(outrot != NULL){
		fclose(outrot);
	}

	return;
}
