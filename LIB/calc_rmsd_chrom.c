#include "gaboom.h"

/******************************************************************************
 * SUBROUTINE calc_rmsd_chrom calculates the rmsd between any two chromosomes  
 * present in the population.
 *****************************************************************************/
float calc_rmsd_chrom(FA_Global* FA, GB_Global* GB, const chromosome* chrom, const genlim* gene_lim,atom* atoms,resid* residue,gridpoint* cleftgrid,int npar, int chrom_a, int chrom_b,
                      float* coor_a_dest, float* coor_b_dest, bool calc_rmsd){

	float rmsd_chrom=0.0f;
	int i = 0,j = 0,k = 0,l = 0,m = 0;
	int cat;
	int rot;
    
	//float* coor_a;
	//float* coor_b;
    
    float coor_a[MAX_ATM_HET*3];
    float coor_b[MAX_ATM_HET*3];
    
    if(coor_a_dest == NULL){
        coor_a_dest = coor_a;
    }
    
    if(coor_b_dest == NULL){
        coor_b_dest = coor_b;
    }
    
	uint grd_idx;
	int normalmode=-1;
	int rot_idx=0;

	for(k=0;k<2;k++)
	{
		j=chrom_a;
		if(k==1){j=chrom_b;}
        
		for(i=0;i<npar;i++){ FA->opt_par[i] = chrom[j].genes[i].to_ic; }
  
		/*
		  printf("%2d (",j);
		  for(l=0;l<GB->num_genes;l++) printf("%12.6f ",FA->opt_par[l]);
		  printf(") ");
		  printf(" value=%11.6f fitnes=%11.6f\n",chrom[j].evalue,chrom[j].fitnes);
		*/
    
    
		for(i=0;i<npar;i++)
		{
			//printf("[%8.3f]",FA->opt_par[i]);
      
			if(FA->map_par[i].typ==-1) 
			{ //by index
				grd_idx = (uint)FA->opt_par[i];
				//printf("FA->opt_par(index): %d\n", grd_idx);
				//PAUSE;
				atoms[FA->map_par[i].atm].dis = cleftgrid[grd_idx].dis;
				atoms[FA->map_par[i].atm].ang = cleftgrid[grd_idx].ang;
				atoms[FA->map_par[i].atm].dih = cleftgrid[grd_idx].dih;
	
			}
			else if(FA->map_par[i].typ==0) 
			{
				atoms[FA->map_par[i].atm].dis = (float)FA->opt_par[i];
	
			}
			else if(FA->map_par[i].typ==1) 
			{
				atoms[FA->map_par[i].atm].ang = (float)FA->opt_par[i];
			}
			else if(FA->map_par[i].typ==2)
			{
				atoms[FA->map_par[i].atm].dih = (float)FA->opt_par[i];
	
				j=FA->map_par[i].atm;
				cat=atoms[j].rec[3];
				if(cat != 0){
					while(cat != FA->map_par[i].atm){
						atoms[cat].dih=atoms[j].dih + atoms[cat].shift; 
						j=cat;
						cat=atoms[j].rec[3];
					}
				}
			}else if(FA->map_par[i].typ==3)
			{ //by index
				grd_idx = (uint)FA->opt_par[i];
				//printf("FA->opt_par(index): %d\n", grd_idx);
				//PAUSE;
	
				// serves as flag , but also as grid index
				normalmode=grd_idx;
	
			}else if(FA->map_par[i].typ==4)\
			{
				rot_idx = (int)(FA->opt_par[i]+0.5);
	
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
  
		/* rebuild cartesian coordinates of optimized residues*/
		for(i=0;i<FA->nors;i++){
			/*printf("nors=%d opt_res[%d]=%d nmov[%d]=%d\n",
			  i,i,FA->opt_res[i],i,FA->nmov[i]);
			  for(j=0;j<FA->nmov[i];j++){
			  printf("mov[%d][%d]=%d\n",i,j,FA->mov[i][j]);
			  }
			  PAUSE;*/
			buildcc(FA,atoms,FA->nmov[i],FA->mov[i]);
		}

		// residue that is optimized geometrically (ligand)
		l=atoms[FA->map_par[0].atm].ofres;

		m=0;
		rot=residue[l].rot;
		for(i=residue[l].fatm[rot];i<=residue[l].latm[rot];i++){
			//printf("i:%d %f %f %f\n",i,atoms[i].coor[0],
			//     atoms[i].coor[1],
			//     atoms[i].coor[2]);
			for(j=0;j<3;j++)
			{
				if(k==0) coor_a_dest[m*3+j]=atoms[i].coor[j];
				if(k==1) coor_b_dest[m*3+j]=atoms[i].coor[j];
			}
			m++;
		}    
	}
  
    if(calc_rmsd){
        for(i=0;i<m;i++)
            rmsd_chrom += sqrdist(&coor_a_dest[i*3],&coor_b_dest[i*3]);
        
        rmsd_chrom = sqrt(rmsd_chrom/((float)m));
    }
	//printf("RMSD=%f\n",rmsd_chrom);
	//PAUSE;

	return rmsd_chrom;
}
