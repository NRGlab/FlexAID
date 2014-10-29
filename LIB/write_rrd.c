#include "gaboom.h"
#include "boinc.h"

/******************************************************************************
 * SUBROUTINE write_pdb writes a PDB file
 ******************************************************************************/
int write_rrd(FA_Global* FA,GB_Global* GB,const chromosome* chrom, const genlim* gene_lim, atom* atoms,resid* residue,gridpoint* cleftgrid,int* Clus_GAPOP,float* Clus_RMSDT,char outfile[]){
	FILE *outfile_ptr;
	int  i,j,l;
	char sufix[10];
	char tmp_end_strfile[MAX_PATH__];
	float rmsd;

	sprintf(sufix,".rrd");
	strcpy(tmp_end_strfile,outfile);
	strcat(tmp_end_strfile,sufix);
	
	outfile_ptr=NULL;
	if(!OpenFile_B(tmp_end_strfile,"w",&outfile_ptr)){
		Terminate(6);
	}else{
		for(j=0;j<GB->num_chrom;j++){
			for(i=0;i<GB->num_genes;i++){
				FA->opt_par[i] = chrom[j].genes[i].to_ic;
			}
			
			rmsd=calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par);
			// rmsd=calc_Hungarian_RMSD(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par);
			fprintf(outfile_ptr,"%3d %3d %8.5f %8.5f %8.5f [",j,Clus_GAPOP[j],Clus_RMSDT[j],rmsd,chrom[j].evalue);
			for(l=0;l<FA->npar;l++){fprintf(outfile_ptr,"%8.5f ",FA->opt_par[l]);}
			fprintf(outfile_ptr,"]\n");
		}
	}

	CloseFile_B(&outfile_ptr,"w");
	return(0);
}
