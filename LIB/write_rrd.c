#include "gaboom.h"
#include "boinc.h"

/******************************************************************************
 * SUBROUTINE write_pdb writes a PDB file
 ******************************************************************************/
int write_rrd(FA_Global* FA,GB_Global* GB,const chromosome* chrom, const genlim* gene_lim, atom* atoms,resid* residue,gridpoint* cleftgrid, int* Clus_GAPOP,float* Clus_RMSDT,char outfile[]){
	FILE *outfile_ptr;
	int  i,j,l;
	char sufix[10];
	char tmp_end_strfile[MAX_PATH__];
	float rmsd = 0.0f;
	float rmsd_corrected = 0.0f;
	bool Hungarian = false;

	sprintf(sufix,".rrd");
	strcpy(tmp_end_strfile,outfile);
	strcat(tmp_end_strfile,sufix);
	
	outfile_ptr=NULL;
	if(!OpenFile_B(tmp_end_strfile,"w",&outfile_ptr)){
		Terminate(6);
	}else{
		for(j=0;j<GB->num_chrom;j++)
		{
			for(i=0;i<GB->num_genes;i++)
			{
				FA->opt_par[i] = chrom[j].genes[i].to_ic;
			}

			// Calculating standard RMSD saved inton 'rmsd'
			Hungarian = false;
			rmsd = calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian);

			// Calculating symmetry-corrected RMSD (with the help of the Hungarian) saved into 'rmsd_corrected'
			Hungarian = true;
			rmsd_corrected = calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian);
			
			fprintf(outfile_ptr,"%3d %3d %8.5f %8.5f %8.5f %8.5f [ ",j,Clus_GAPOP[j],Clus_RMSDT[j],rmsd,rmsd_corrected,chrom[j].evalue);
			for(l=0;l<FA->npar;l++){fprintf(outfile_ptr,"%8.5f ",FA->opt_par[l]);}
			fprintf(outfile_ptr," ]\n");
		}
	}

	CloseFile_B(&outfile_ptr,"w");
	return(0);
}

int write_DensityPeak_rrd(FA_Global* FA, GB_Global* GB, const chromosome* chrom, const genlim* gene_lim, atom* atoms, resid* residue, gridpoint* cleftgrid, ClusterChrom* Chrom, DPcluster* Clust, float* RMSD, char outfile[])
{
	FILE *outfile_ptr;
	int i,j,k,l;
	char sufix[10];
	char tmp_end_strfile[MAX_PATH__];
	float rmsd = 0.0f;
	float rmsd_corrected = 0.0f;
	float ClusRMSD = 0.0f;
	bool Hungarian = false;

	sprintf(sufix,".rrd");
	strcpy(tmp_end_strfile, outfile);
	strcat(tmp_end_strfile, sufix);

	outfile_ptr = NULL;
	if(!OpenFile_B(tmp_end_strfile,"w",&outfile_ptr)) 
	{
		Terminate(6);
	}
	else
	{
		for(j=0; j<GB->num_chrom;++j)
		{ 
			for(i=0;i<GB->num_genes;++i)
			{
				FA->opt_par[i] = chrom[j].genes[i].to_ic;
			}

			Hungarian = false;
			rmsd = calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian);

			// Calculating symmetry-corrected RMSD (with the help of the Hungarian) saved into 'rmsd_corrected'
			Hungarian = true;
			rmsd_corrected = calc_rmsd(FA,atoms,residue,cleftgrid,FA->npar,FA->opt_par, Hungarian);
			for(k=j+1;k<GB->num_chrom;++k)
			{
				if(Chrom[k].Cluster == Chrom[j].Cluster && Chrom[k].isCenter == true)
				{
					ClusRMSD = RMSD[K(Chrom[k].index, Chrom[j].index, GB->num_chrom)];
				}
			}
			fprintf(outfile_ptr, "%3d %3d %8.5f %8.5f %8.5f %8.5f [", j, Chrom[j].Cluster, ClusRMSD, rmsd, rmsd_corrected, chrom[j].evalue);
			for(l=0; l<FA->npar; ++l) fprintf(outfile_ptr,"%8.5f ",FA->opt_par[l]);
			fprintf(outfile_ptr, "]\n" );
		}
	}

	CloseFile_B(&outfile_ptr,"w");
	return(0);
}