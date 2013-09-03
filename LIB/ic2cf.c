#include "gaboom.h"
#include "boinc.h"

/******************************************************************************
 * SUBROUTINE ic2cf gets a vector with internal coordinates rebuilds the 
 * cartesian coordinates and calculates the complementarity function. Its 
 * input vector has the list of ic's that are to be optimized and a global
 * vector contains the information of what kind of variable each item in icv
 * is and to which residue it belongs.
 *****************************************************************************/

//THE PROCEDURE SHOULD RECEIVE A 2ND SET OF GENES THAT ENCODES FOR THE ROTAMER DISTRIBUTION IN THE BPK

cfstr ic2cf(FA_Global* FA,VC_Global* VC,atom* atoms,resid* residue,gridpoint* cleftgrid,int npar, double* icv){
  
	static int nbranch = 0;
	
	int i,j,k;
	int cat;    /* atom number constrained to the one considered */

	cfstr cf;
	static cfstr cf_clash = { 0.0, 0.0, CLASH_PENALTY_VALUE, 0.0, 1 };
	
	int rclash=0;

	psFlexDEE_Node psFlexDEENode;
	int dee_val;
	int rotflag;
  
	unsigned int grd_idx;
	unsigned int rot_idx;
	//float threshold=10.0;
	//int nflxchk=-1;
	int normalmode=-1;
	int deelig_list[100];
	
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
	if(normalmode > -1){
		alter_mode(atoms,residue,FA->normal_grid[normalmode],FA->res_cnt,FA->normal_modes);
	}
	
	/* rebuild cartesian coordinates of optimized residues*/
  
  
	for(i=0;i<FA->nors;i++){ //number of optimized residues
		//printf("nors=%d opt_res[%d]=%d nmov[%d]=%d\n",i,i,FA->opt_res[i],i,FA->nmov[i]);
		//for(j=0;j<nmov[i];j++){printf("mov[%d][%d]=%d\n",i,j,FA->mov[i][j]);}
		//PAUSE;
    
		buildcc(FA,atoms,FA->nmov[i],FA->mov[i]);
	}
  
	vector< pair<int,int> > intraclashes;
	if(vcfunction(FA,VC,atoms,residue,intraclashes)){
		return cf_clash;
	}
	
	cf.com = 0.0;
	cf.wal = 0.0;
	cf.sas = 0.0;
	cf.con = 0.0;
	cf.rclash = 0;
    
	for(i=0;i<FA->num_optres;i++){
    
		resid* res = &residue[FA->optres[i].rnum];
		
		// flexible side-chain optimization
		if ( !FA->optres[i].type ) {
  
			if ( FA->optres[i].cf.rclash == 1 ) { 

				/*
				  printf("%s %c %d is clashing\n",
				  residue[FA->optres[i].rnum].name,
				  residue[FA->optres[i].rnum].chn,
				  residue[FA->optres[i].rnum].number);
				*/
	
				rclash = 1; 
				
			}
      
		}else{
			
			//int fatm = res->fatm[0];
			
			vector< pair<int,int> >::iterator it;
			for(it=intraclashes.begin(); it!=intraclashes.end(); ++it)
			{
				for(k=1; k<=res->fdih; k++){
					deelig_list[k] = -1000;
				}
				
				// flex bonds list
				int fbindex = 0;
				int* fblist = res->shortflex[it->first][it->second];
				
				//printf("between[%d][%d]\n", atoms[it->first+fatm].number, atoms[it->second+fatm].number);
				//cout << fblist[fbindex] << endl;
				
				while(fblist[fbindex] != -1){
					if(atoms[res->bond[fblist[fbindex]]].par != NULL){
						deelig_list[fblist[fbindex]] =
							(int)(atoms[res->bond[fblist[fbindex]]].dih + 0.5);
					}
					fbindex++;
				}
				
				struct deelig_node_struct* node = FA->deelig_root_node;
				
				bool add = false;
				for(k=1; k<=res->fdih; k++){
					std::map<int, struct deelig_node_struct*>::iterator it;
					it = node->childs.find(deelig_list[k]);
					
					if(it == node->childs.end()){
						struct deelig_node_struct* deelig_child_node = new struct deelig_node_struct;
						if(!deelig_child_node){
							fprintf(stderr, "ERROR: memory allocation error for deelig_child_node\n");
							Terminate(2);
						}
						
						//if(k==1) cout << "new node added " << deelig_list[k] << endl;
						node->childs[deelig_list[k]] = deelig_child_node;
						
						deelig_child_node->parent = node;
						node = deelig_child_node;
						add = true;
					}else{
						node = it->second;
					}
				}
				
				if(add) nbranch++;
				//cout << "total branches " << nbranch << endl;
			}
			
		}
    
		/*
		  printf("optres[%2d].cf  .wal = %.3f\n               .com = %.3f\n               .sas = %.3f\n               .con = %.3f\n",
		  i,FA->optres[i].cf.wal,FA->optres[i].cf.com,FA->optres[i].cf.sas,FA->optres[i].cf.con);
		*/
        
		//sum += (FA->optres[i].cf.com - FA->optres[i].cf.wal + FA->optres[i].cf.sas - FA->optres[i].cf.con);
    
		cf.com += FA->optres[i].cf.com;
		cf.wal += FA->optres[i].cf.wal;
		cf.sas += FA->optres[i].cf.sas;
		cf.con += FA->optres[i].cf.con;

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

 
	return cf;

}

void add_intraclash_key(void)
{

}


#ifdef _WIN32
double get_apparent_cf_evalue(cfstr* cf) {
#else
	double get_apparent_cf_evalue(cfstr* cf) {
#endif
		return cf->com + cf->wal + cf->sas;
	}
    
#ifdef _WIN32
	double get_cf_evalue(cfstr* cf) {
#else
		double get_cf_evalue(cfstr* cf) {
#endif
			return cf->com + cf->wal + cf->sas + cf->con;
		}
