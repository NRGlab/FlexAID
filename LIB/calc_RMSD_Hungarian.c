#include "flexaid.h"

/********************************************************************************
 * This function calculates the RSMD between atomic coordinates of the atoms in *
 * the register ori_ligatm and those for the atoms of the ligand in             *
 * residue[opt_res[0]] after reconstructing the coordinates using opt_par       *
 ********************************************************************************/

float calc_Hungarian_RMSD(FA_Global* FA, atom* atoms, resid* residue, gridpoint* cleftgrid,int npar, const double* icv)
{
	int i,j,k;
 	float rmsd = 0.0f;
 	
	int nTypes = 0;
	int nUniqueTypes = 0;
 	bool unique_flag = false;

 	// Counting the number of unique atom types
 	// Initilization of the 'unique_atom_type[]'
 	int unique_atom_type[FA->num_het_atm];
 	for(i=0; i<FA->num_het_atm;i++)
 	{
 		unique_atom_type[i] = 0;
 	}

	for(i = 0; i<FA->num_het; i++)
	{
		for(j = 1; j<=FA->res_cnt; j++)
		{
			if(residue[j].number == FA->het_res[i])
			{
				// The 'for loop' below iterates through all the atoms of the ligand
 				for(k=residue[j].fatm[0]; k<=residue[j].latm[0]; k++)
 				{
					for(int jj=0; jj < nUniqueTypes+1; jj++)
					{
						if(atoms[k].type == unique_atom_type[jj])
						{
							unique_flag = false;
							break;
						}
					}
					if(unique_flag)
					{
						unique_atom_type[nUniqueTypes] = atoms[k].type;
						nUniqueTypes++;
					}
					unique_flag = true;
 				}
 				// Counting the number of unique atom types
	 			for(int z = 0; z < nUniqueTypes; z++)
	 			{
	 				for(k=residue[j].fatm[0]; k<=residue[j].latm[0]; k++)
	 				{
	 					if(atoms[k].type == unique_atom_type[z])
	 					{
	 						nTypes++;
	 					}
	 				}
	 			
	 				// Initializing Hungarian algorithm variables
					float matrix[nTypes][nTypes];
					float matrix_original[nTypes][nTypes];
					int matrix_cases[nTypes][nTypes];
					
					int row_count[nTypes];
					int column_count[nTypes];
					int column_assigned[nTypes];
					int matrix_match[nTypes];
	 				
	 				int count_i = -1;
					int count_j = -1;

					// Step 1 -> Reducing the matrix
					for(k=residue[j].fatm[0]; k<=residue[j].latm[0]; k++) // Foreach Atom in Lig
					{
						if(atoms[k].type == unique_atom_type[z])
						{
							count_i++;
							count_j = -1;
							for(int l = 0; l < FA->num_het_atm; l++)
							{
								if(atoms[l].type == unique_atom_type[z])
								{
									
								}
							}
						}
					}
				}
			}
		}
	}
}