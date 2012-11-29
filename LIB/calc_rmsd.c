#include "flexaid.h"

/********************************************************************************
 * This function calculates the RSMD between atomic coordinates of the atoms in *
 * the register ori_ligatm and those for the atoms of the ligand in             *
 * residue[opt_res[0]] after reconstructing the coordinates using opt_par       *
 ********************************************************************************/

float calc_rmsd(FA_Global* FA,atom* atoms,resid* residue, gridpoint* cleftgrid,int npar, const double* icv){
	float rmsd=0.0f;
	int i,cat,j,l;
	uint grd_idx;
	int normalmode=-1;
	int rot_idx=0,rot;
  
	for(i=0;i<npar;i++){
		//printf("[%8.3f]",icv[i]);
    
		if(FA->map_par[i].typ==-1) { //by index
			grd_idx = (uint)icv[i];
			//printf("icv(index): %d\n", grd_idx);
			//PAUSE;
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
			rot_idx = (int)(icv[i]+0.5);
      
			residue[atoms[FA->map_par[i].atm].ofres].rot=rot_idx;
      
			/*
			  printf("residue[%d].rot[%d] - fatm=%d - latm=%d\n",
			  residue[atoms[FA->map_par[i].atm].ofres].number,
			  residue[atoms[FA->map_par[i].atm].ofres].rot,
			  residue[atoms[FA->map_par[i].atm].ofres].fatm[rot_idx],
			  residue[atoms[FA->map_par[i].atm].ofres].latm[rot_idx]);
			*/
      
		}
    
	}

  
	if(normalmode > -1)
		alter_mode(atoms,residue,FA->normal_grid[normalmode],FA->res_cnt,FA->normal_modes);
  
	// rebuild cartesian coordinates of optimized residues
	for(i=0;i<FA->nors;i++){
		/*printf("nors=%d opt_res[%d]=%d nmov[%d]=%d\n",
		  i,i,FA->opt_res[i],i,FA->nmov[i]);
		  for(j=0;j<FA->nmov[i];j++){
		  printf("mov[%d][%d]=%d\n",i,j,FA->mov[i][j]);
		  }
		  PAUSE;
		*/
		buildcc(FA,atoms,FA->nmov[i],FA->mov[i]);
	}
 
	l=0;
	for(i=1; i<=FA->res_cnt; i++){
		rot = residue[i].rot;
		for(j=residue[i].fatm[rot]; j<=residue[i].latm[rot]; j++){
			if(atoms[j].coor_ref != NULL){
				rmsd += sqrdist(atoms[j].coor,atoms[j].coor_ref);
				l++;
			}
		}
	}
	
	rmsd = sqrtf(rmsd/((float)l));
	
	//printf("RMSD=%f\n",rmsd); 
	//PAUSE;
 
	return rmsd;
}
