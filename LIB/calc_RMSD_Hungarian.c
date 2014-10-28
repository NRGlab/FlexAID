#include "flexaid.h"

float calc_Hungarian_RMSD(FA_Global* FA, atom* atoms, resid* residue, gridpoint* cleftgrid,int npar, const double* icv)
{
	int i,j,k;
 	float rmsd = 0.0f;
 	float total_assignment = 0;
	int nTypes = 0;
	int nUniqueTypes = 0;
 	bool unique_flag = false;

 	// Counting the number of unique atom types
 	// Initilization of the 'unique_atom_type[]'
 	int *unique_atom_type;
 	unique_atom_type = (int*) malloc(sizeof(int) * FA->num_het_atm);
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
	 			
	 				// Declaration and memomry allocation for the Hungarian algorihm arrays 
					// ~initialize();
					int **matrix_case;
					float **matrix, **matrix_original;
					matrix = (float**) malloc(sizeof(float *) * nTypes);
					matrix_original = (float**) malloc(sizeof(float *) * nTypes);
					matrix_case = (int**) malloc(sizeof(int*) * nTypes);
					
					for(k = 0; k < nTypes; k++)
					{
						matrix[k] = (float *) malloc(sizeof(float) * nTypes);
						matrix_original[k] = (float *) malloc(sizeof(float) * nTypes);
						matrix_case[k] = (int *) malloc(sizeof(int) * nTypes);
					}
					
					int *row_count, *column_count, *row_assigned, *column_assigned, *matrix_match;
					row_count = (int*) malloc(sizeof(int) * nTypes);
					column_count = (int*) malloc(sizeof(int) * nTypes);
					row_assigned = (int*) malloc(sizeof(int) * nTypes);
					column_assigned = (int*) malloc(sizeof(int) * nTypes);
					matrix_match = (int*) malloc(sizeof(int) * nTypes);
	 				
	 				int count_i = -1;
					int count_j = -1;

					// Filling the matrix
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
									count_j++;
									matrix[count_i][count_j] = sqrdist(atoms[k].coor,atoms[l].coor_ref);
								}
							}
						}
					}
					Hungarian(matrix, matrix_original, matrix_case, unique_atom_type, row_count, column_count, row_assigned, column_assigned, matrix_match, nTypes, FA->num_het_atm);
					float assignment = 0;
					for(int ii = 0; ii < nTypes; ii++)
					{
						assignment += matrix_original[ii][matrix_match[ii]];
					}
					total_assignment += assignment;
					
					// Free dynamically allocated memory
					if(matrix != NULL)
					{
						for(int ii = 0; ii < nTypes; ii++) {free(matrix[ii]);}
						free(matrix);
					}
					if(matrix_original != NULL)
					{
						for(int ii = 0; ii < nTypes; ii++) {free(matrix_original[ii]);}
						free(matrix_original);
					}
					if(matrix_case != NULL)
					{
						for(int ii = 0; ii < nTypes; ii++) {free(matrix_case[ii]);}
						free(matrix_case);
					}
					if(row_count != NULL) {free(row_count);}
					if(column_count != NULL) {free(column_count);}
					if(row_assigned != NULL) {free(row_assigned);}
					if(column_assigned != NULL) {free(column_assigned);}
					if(matrix_match != NULL) {free(matrix_match);}
				}
				rmsd = sqrt(total_assignment/FA->num_het_atm);

				// Freeing the unique_atom_type[] dynamically allocated memory (check if it could be placed elsewhere)
				if(unique_atom_type != NULL) {free(unique_atom_type);}
				
				// Return symetry corrected RMSD
				return(rmsd);
			}
		}
	}
}

void Hungarian(float** matrix, float** matrix_original, int** matrix_case, int* unique_atom_type, int* row_count, int* column_count, int* row_assigned, int* column_assigned, int* matrix_match, int nTypes, int num_het_atoms)
{
	// Step 1. Copy and Reduce the matrix
	Hungarian_reduce_matrix(matrix, matrix_original, nTypes);

	int cycle_counter = 0;
	while( cycle_counter < (num_het_atoms * num_het_atoms) )
	{
		// Step 2. Attemp to assign one job to each worker
		int number_assigned = Hungarian_assign_jobs(matrix, matrix_original, matrix_case, unique_atom_type, row_count, column_count, row_assigned, column_assigned, matrix_match, nTypes, num_het_atoms);
		if(number_assigned == nTypes) {break;}
		// Step 3. Draw the minimum number of lines required to cross out all 0s from matrix[][]
		Hungarian_draw_line(matrix, matrix_original, matrix_case, row_count, column_count, row_assigned, column_assigned, matrix_match, nTypes);
		// Step 4. Update cells of matrix[][] from lines drawn
		Hungarian_update_matrix(matrix, matrix_case, nTypes);
		cycle_counter++;
	}
}
int Hungarian_assign_jobs(float** matrix, float** matrix_original, int** matrix_case, int* unique_atom_type, int* row_count, int* column_count, int* row_assigned, int* column_assigned, int* matrix_match, int nTypes, int num_het_atoms)
{
	Hungarian_reset_match(matrix_match, nTypes);
	Hungarian_reset_assigned_row(row_assigned, nTypes);
	Hungarian_reset_assigned_column(column_assigned, nTypes);
	int number_assigned = 0;
	while(number_assigned < nTypes+1)
	{
		Hungarian_reset_row_count(row_count, nTypes);
		Hungarian_reset_column_count(column_count, nTypes);

		// 1. Count the number of 0s in each row and column
		for(int i=0; i<nTypes; i++)
		{
			for(int j=0; j<nTypes;j++)
			{
				if(matrix[i][j] == 0)
				{
					if(row_assigned[i] == 0 && column_assigned[j] == 0)
					{
						row_count[i]++;
						column_count[j]++;
					}
				}
			}
		}

	}
}
void Hungarian_reset_match(int* matrix_match, int nTypes)
{

}
void Hungarian_reset_assigned_row(int* row_assigned, int nTypes)
{

}
void Hungarian_reset_assigned_column(int* column_assigned, int nTypes)
{

}
void Hungarian_reset_row_count(int* row_count, int nTypes)
{

}
void Hungarian_reset_column_count(int* column_count, int nTypes)
{

}
void Hungarian_update_matrix(float** matrix, int** matrix_case, int nTypes)
{

}
void Hungarian_draw_line(float** matrix, float** matrix_original, int** matrix_case, int* row_count, int* column_count, int* row_assigned, int* column_assigned, int* matrix_match, int nTypes)
{

}
void Hungarian_reduce_matrix(float** matrix, float** matrix_original, int nTypes)
{
	for(int ii=0; ii<nTypes; ii++)
	{
		for(int jj=0; jj < nTypes; jj++)
		{
			matrix_original[ii][jj] = matrix[ii][jj];
		}
	}

	for(int ii=0; ii<nTypes; ii++)
	{
		float minimal_value = float.MaxValue;
		for(int jj=0; jj < nTypes; jj++)
		{
			if(matrix[ii][jj] < minimal_value)
			{
				minimal_value = matrix[ii][jj];
			}
		}
		for(int jj =0; jj < nTypes; jj++)
		{
			matrix[ii][jj] -= minimal_value;
		}
	}

	for(int jj=0; jj<nTypes; jj++)
	{
		float minimal_value = float.MaxValue;
		for(int ii =0; ii<nTypes; ii++)
		{
			if(matrix[ii][jj] < minimal_value)
			{
				minimal_value = matrix[ii][jj];
			}
		}
		for(int ii=0; ii<nTypes;ii++)
		{
			matrix[ii][jj] -= minimal_value;
		}
	}
	return;
}