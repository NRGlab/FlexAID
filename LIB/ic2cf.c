#include "gaboom.h"

/******************************************************************************
 * SUBROUTINE ic2cf gets a vector with internal coordinates rebuilds the 
 * cartesian coordinates and calculates the complementarity function. Its 
 * input vector has the list of ic's that are to be optimized and a global
 * vector contains the information of what kind of variable each item in icv
 * is and to which residue it belongs.
 *****************************************************************************/

//THE PROCEDURE SHOULD RECEIVE A 2ND SET OF GENES THAT ENCODES FOR THE ROTAMER DISTRIBUTION IN THE BPK

double ic2cf(FA_Global* FA,VC_Global* VC,atom* atoms,resid* residue,gridpoint* cleftgrid,int npar, double* icv){
  
	int i,j,k;
	int cat;    /* atom number constrained to the one considered */

	double sum=0.0;
	int rclash=0;

	psFlexDEE_Node psFlexDEENode;
	int dee_val;
	int rotflag;
  
	unsigned int grd_idx;
	unsigned int rot_idx;
	//float threshold=10.0;
	//int nflxchk=-1;
	int normalmode=-1;

	//int rigid_clash=0;
	//float rand=0.0;
	//float min_dis=10.0;
  
	// copy values from icv into respective srtructure atom ic fields 
	// andcompute the ic of a constrained atom prior to reconstruction

	//for(i=0;i<npar;i++){printf("[%8.3f]",icv[i]);}printf("\n");
	//PAUSE;

	//printf("NEW INDIVIDUAL=");
	for(i=0;i<npar;i++){
		//printf("[%8.3f]",icv[i]);

		if(FA->map_par[i].typ==-1) { //by index

			grd_idx = (uint)icv[i];
			atoms[FA->map_par[i].atm].dis = cleftgrid[grd_idx].dis;
			atoms[FA->map_par[i].atm].ang = cleftgrid[grd_idx].ang;
			atoms[FA->map_par[i].atm].dih = cleftgrid[grd_idx].dih;
      
		}else if(FA->map_par[i].typ==0)  {
			atoms[FA->map_par[i].atm].dis = (float)icv[i];
      
		}else if(FA->map_par[i].typ==1)  {
			atoms[FA->map_par[i].atm].ang = (float)icv[i];
      
		}else if(FA->map_par[i].typ==2)  {
			atoms[FA->map_par[i].atm].dih = (float)icv[i];
      
			j=FA->map_par[i].atm;
			cat=atoms[j].rec[3];
			if(cat != 0){
				while(cat != FA->map_par[i].atm){
					atoms[cat].dih=atoms[j].dih + atoms[cat].shift; 
					j=cat;
					cat=atoms[j].rec[3];
				}
			}
            
		}else if(FA->map_par[i].typ==3) { //by index
			grd_idx = (uint)icv[i];
			//printf("icv(index): %d\n", grd_idx);
			//PAUSE;
      
			// serves as flag , but also as grid index
			normalmode=grd_idx;
      
		}else if(FA->map_par[i].typ==4)  {
			rot_idx = (uint)(icv[i]+0.5f);
      
			residue[atoms[FA->map_par[i].atm].ofres].rot=(int)rot_idx;
      
			/*
			  printf("residue[%d].rot[%d] - fatm=%d - latm=%d\n",
			  residue[atoms[FA->map_par[i].atm].ofres].number,
			  residue[atoms[FA->map_par[i].atm].ofres].rot,
			  residue[atoms[FA->map_par[i].atm].ofres].fatm[rot_idx],
			  residue[atoms[FA->map_par[i].atm].ofres].latm[rot_idx]);
			*/
      
		}
    
	}
	//printf("HERE\n");
	//PAUSE;
  
	// do not alter default (ini) protein conf.
	if(normalmode > -1)
		alter_mode(atoms,residue,FA->normal_grid[normalmode],FA->res_cnt,FA->normal_modes);
  
	/* rebuild cartesian coordinates of optimized residues*/
  
  
	for(i=0;i<FA->nors;i++){ //number of optimized residues
		//printf("nors=%d opt_res[%d]=%d nmov[%d]=%d\n",i,i,FA->opt_res[i],i,FA->nmov[i]);
		//for(j=0;j<nmov[i];j++){printf("mov[%d][%d]=%d\n",i,j,FA->mov[i][j]);}
		//PAUSE;
    
		buildcc(FA,atoms,FA->nmov[i],FA->mov[i]);
	}
  
  
	if (strcmp(FA->complf,"VCT")==0){
		if(vcfunction(FA,VC,atoms,residue)){
			return CLASH_PENALTY_VALUE;
		}
	}else if (strcmp(FA->complf,"SPH")==0){
		if(spfunction(FA,atoms,residue)){
			return CLASH_PENALTY_VALUE;
		}
	}
  

	for(i=0;i<FA->num_optres;i++){
    
		// flexible side-chain optimization
		if ( i != 0 ) {
  
			if ( FA->optres[i].cf.rclash == 1 ) { 

				/*
				  printf("%s %c %d is clashing\n",
				  residue[FA->optres[i].rnum].name,
				  residue[FA->optres[i].rnum].chn,
				  residue[FA->optres[i].rnum].number);
				*/
	
				rclash = 1; 

			}
      
		}
    
		/*
		  printf("optres[%2d].cf  .wal = %.3f\n               .com = %.3f\n               .sas = %.3f\n               .con = %.3f\n",
		  i,FA->optres[i].cf.wal,FA->optres[i].cf.com,FA->optres[i].cf.sas,FA->optres[i].cf.con);
		*/

		sum += (FA->optres[i].cf.com - FA->optres[i].cf.wal + FA->optres[i].cf.sas - FA->optres[i].cf.con);
    
	}

  
	// add rotamer list to dee list
	if (FA->useflexdee > 0 && rclash) {
    
		NEW( psFlexDEENode, sFlexDEE_Node );

		psFlexDEENode->rotlist = (int*)malloc(FA->nflxsc_real*sizeof(int));
    
    
		// fill rotamer list
		k=0;
		rotflag=0;
    
		for(j=0;j<FA->nflxsc;j++){
      
			if(residue[FA->flex_res[j].inum].trot > 0   &&
			   FA->flex_res[j].cflag != 0){
	
				psFlexDEENode->rotlist[k++] = residue[FA->flex_res[j].inum].rot;
	
				if(residue[FA->flex_res[j].inum].rot != 0) { rotflag=1; }
			}
      
		}
    
    
		// do not add initial conformation to DEE list
		if ( rotflag ) {
      
			/*
			  printf("\n-----------------\nCreating new node...\n");
			  printf("DEE to add = ");for(k=0;k<FA->nflxsc_real;k++){printf("%3d",psFlexDEENode->rotlist[k]);}printf("\n");
	
			  dee_print(FA->psFlexDEENode,FA->nflxsc_real);
	
			  //getchar();
			  */
      
    
			if( FA->psFlexDEENode ) {
	
				//FA->psFlexDEENode = FA->psFlexDEENode->last;
	
				while ( FA->psFlexDEENode->next != NULL ) {
					FA->psFlexDEENode = FA->psFlexDEENode->next;
				}
	
				dee_val = dee_pivot(psFlexDEENode,&FA->psFlexDEENode,1,FA->FlexDEE_Nodes,(int)((FA->FlexDEE_Nodes+1)/2),FA->FlexDEE_Nodes,FA->nflxsc_real);
	
				if ( dee_val == 1 ) {
	  
					if ( FA->psFlexDEENode->next == NULL ) {
	    
						psFlexDEENode->next = NULL;
						psFlexDEENode->prev = FA->psFlexDEENode;
						FA->psFlexDEENode->next = psFlexDEENode;
	    
						psFlexDEENode->first = FA->psFlexDEENode;
	    
						dee_last(FA->psFlexDEENode,psFlexDEENode);
	    
					} else {
	    
						psFlexDEENode->first = FA->psFlexDEENode->first;
						psFlexDEENode->last = FA->psFlexDEENode->last;
	    
						psFlexDEENode->next = FA->psFlexDEENode->next;
						psFlexDEENode->prev = FA->psFlexDEENode;
						FA->psFlexDEENode->next = psFlexDEENode; 
						psFlexDEENode->next->prev = psFlexDEENode;
	    
					}
	  
					FA->FlexDEE_Nodes++;
	  
				} else if ( dee_val == -1 ) {
	  
					if ( FA->psFlexDEENode->prev == NULL ) {
	    
						psFlexDEENode->prev = NULL;
						psFlexDEENode->next = FA->psFlexDEENode;
						FA->psFlexDEENode->prev = psFlexDEENode;
	    
						psFlexDEENode->last = FA->psFlexDEENode;
	    
						dee_first(FA->psFlexDEENode,psFlexDEENode);
	    
					} else {
	    
						psFlexDEENode->first = FA->psFlexDEENode->first;
						psFlexDEENode->last = FA->psFlexDEENode->last;
	    
						psFlexDEENode->prev = FA->psFlexDEENode->prev;
						psFlexDEENode->next = FA->psFlexDEENode;
						FA->psFlexDEENode->prev = psFlexDEENode; 
						psFlexDEENode->prev->next = psFlexDEENode;
	    
					}
	  
					FA->FlexDEE_Nodes++;
	  
				} else {
	  
					FREE(psFlexDEENode);
	  
				}
	
			} else {
	
				FA->psFlexDEENode = psFlexDEENode;
	
				FA->psFlexDEENode->next = NULL;
				FA->psFlexDEENode->prev = NULL;
	
				FA->psFlexDEENode->first = FA->psFlexDEENode;
				FA->psFlexDEENode->last = FA->psFlexDEENode;
	
				FA->FlexDEE_Nodes++;
	
			}
      
		} 
    
	}

 
	return sum;

}
