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
					
					// FREEING DYNAMICALLY ALLOCATED MEMORY below
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

				// Freeing the unique_atom_type[] dynamically allocated memory (check if it could be placed elsewhere)
				if(unique_atom_type != NULL) {free(unique_atom_type);}
				
				// Return symetry corrected RMSD
				// rmsd = sqrt(total_assignment/FA->num_het_atm);
				// return(rmsd);
			}
		}
	}
	// Return symetry corrected RMSD
	rmsd = sqrt(total_assignment/FA->num_het_atm);
	return(rmsd);
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
	// Reset the matrices to be used in this loop and initilizing number_assigned -> 0
	Hungarian_reset_match(matrix_match, nTypes);
	Hungarian_reset_assigned_row(row_assigned, nTypes);
	Hungarian_reset_assigned_column(column_assigned, nTypes);
	int number_assigned = 0;
	while(number_assigned < nTypes+1)
	{
		// 0. Reset the row_count[] && column_count[]
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

		// 2. Identify the row OR column with the lowest number of 0s
		float minimal_value = 1000000;
		int row_col = 0; // Value of 0 : no more 0s left in matrix; Value of 1 : least number of 0s found in a row; Value of 2 : least number of 0s found in a col 

		int row = -1;
		int col = -1;

		for(int i=0; i<nTypes; i++)
		{
			if(row_count[i] != 0 && row_count[i] < minimal_value)
			{
				minimal_value = row_count[i];
				row = i;
				row_col = 1;
				col = -1;
			}
			if(column_count[i] != 0 && column_count[i] < minimal_value)
			{
				minimal_value = column_count[i];
				col = i;
				row_col = 2;
				row = -1;
			}
		}

		// 3. Make an assignment in the row or column with the least number of 0s
		if(row_col == 0) { break; }
		else if(row_col == 1)
		{
			for(int i=0; i<nTypes; i++)
			{
				if(matrix[row][i] == 0 && column_assigned[i] == 0)
				{
					matrix_match[row] = i;
					row_assigned[row] = 1;
					column_assigned[i] = 1;
					break; // breaks the for loop
				}
			}
		}
		else if(row_col == 2)
		{
			for(int i=0; i<nTypes; i++)
			{
				if(matrix[i][col] == 0 && row_assigned[i] == 0)
				{
					matrix_match[i] = col;
					column_assigned[col] = 1;
					row_assigned[i] = 1;
					break; // breaks the for loop
				}
			}
		}
		number_assigned++;
	} // End of While loop : this will loop until all possible worker / job pairs are formed
	return(number_assigned);
}
void Hungarian_reset_match(int* matrix_match, int nTypes)
{
	for(int i=0; i<nTypes; i++)
	{
		matrix_match[i] = -100;
	}
}
void Hungarian_reset_assigned_row(int* row_assigned, int nTypes)
{
	for(int i=0; i<nTypes; i++)
	{
		row_assigned[i] = 0;
	}
}
void Hungarian_reset_assigned_column(int* column_assigned, int nTypes)
{
	for(int i=0; i<nTypes; i++)
	{
		column_assigned[i] = 0;
	}
}
void Hungarian_reset_row_count(int* row_count, int nTypes)
{
	for(int i=0; i<nTypes; i++)
	{
		row_count[i] = 0;
	}
}
void Hungarian_reset_column_count(int* column_count, int nTypes)
{
	for(int i=0; i<nTypes; i++)
	{
		column_count[i] = 0;
	}
}
void Hungarian_reset_case(int** matrix_case, int nTypes)
{
	for(int i=0;i<nTypes;i++)
	{
		for(int j=0; j<nTypes;j++)
		{
			matrix_case[i][j] = 0;
		}
	}
}
void Hungarian_update_matrix(float** matrix, int** matrix_case, int nTypes)
{
	// 0. Initialization of the variable used to store the minimal value of the matrix
	float minimal_value = 1000000;
	// 1. Find the minimal_value in matrix
	for(int i=0; i<nTypes; i++)
	{
		for(int j=0; j<nTypes; j++)
		{
			if(matrix_case[i][j] == 0)
			{
				if(matrix[i][j] < minimal_value)
				{
					minimal_value = matrix[i][j];
				}
			}
		}
	}
	// 2. Update every element of the matrix[][]
	for(int i=0; i<nTypes; i++)
	{
		for(int j=0; j<nTypes; j++)
		{
			if(matrix_case[i][j] == 0)
			{
				matrix[i][j] -= minimal_value;
			}
			else if(matrix_case[i][j] == 1)
			{
				// Nothing to do here !
			}
			else if(matrix_case[i][j] == 2)
			{
				matrix[i][j] += minimal_value;
			}
		}
	}
	return;
}
void Hungarian_draw_line(float** matrix, float** matrix_original, int** matrix_case, int* row_count, int* column_count, int* row_assigned, int* column_assigned, int* matrix_match, int nTypes)
{
	int lines = 0; // the number of lines required to cross out all 0s in matrix[][]
	// 0. reset the matrices used in this function
	Hungarian_reset_assigned_row(row_assigned, nTypes);
	Hungarian_reset_assigned_column(column_assigned, nTypes);
	Hungarian_reset_case(matrix_case, nTypes);

	// 1. Find/Mark rows without assignment
	for(int i=0; i<nTypes; i++)
	{
		if(matrix_match[i] == -100)
		{
			row_assigned[i] = 1;
		}
	}

	bool flag_update = true; // Flag indicating that an update was made during last iter
	int count_updates = 0;	 // Loop counter

	// Looping until which line to be drawn is decided
	while(flag_update && count_updates < nTypes*2)
	{
		flag_update = false;

		// 2. Mark all columns having 0s in those rows
		for(int i=0; i<nTypes; i++) // Foreach row in matrix
		{
			if(row_assigned[i] == 1) // If this row have been previously assigned
			{
				for(int j=0; j < nTypes; j++) // For each column in this row 
				{
					if(matrix[i][j] == 0) // Check wheter there is a 0
					{
						column_assigned[j] = 1; // Mark && remember this position in colum_assigned
					}
				}
			}
		}

		// 3. Mark all rows having assignments in those coloumns
		for(int i=0; i< nTypes; i++) // for all the column in the matrix[][]
		{
			if(column_assigned[i] == 1) // If a the column have been marked
			{
				for(int j=0; j<nTypes; j++) // Then check each row in this column
				{
					if(matrix_match[j] == i && row_assigned[j] != 1) // and see if an assignment have been made :
					{
						row_assigned[j] = 1; // Mark that row
						flag_update = true;  // Flag the update that have been made
						count_updates++; 	 // Unsure if I need to update the counter here OR elsewhere OR if I need to get rid of it
					} 
				}
			}
		}
	}

	// 4. Draw lines htrough marked columns and unmarked rows
	for(int i=0; i<nTypes; i++)				// for all of the cells in matrix[][]...
	{
		if(row_assigned[i] == 0)			// if a row has not been marked,
		{
			for(int j=0; j<nTypes; j++)		// draw a line through that row by incrementing every entry in
			{
				matrix_case[i][j]++;		// matrix_case[row][] by 1
			}
			lines++;						// a line was just drawn, increment lines by 1
		}
		if(column_assigned[i] == 1)
		{
			for(int j=0; j<nTypes; j++)		// if a column has been marked,
			{								// draw a line through that column by incrementing every entry in
				matrix_case[j][i]++;		// matrix_case[][column] by 1
			}
			lines++;						// a line was just drawn, increment lines by 1
		}
	}
	// every cell in matrix_case[][] will now either be a 0 (no lines covering), 1 (one line covering),
    // or 2 (two lines covering)
	return;
}
void Hungarian_reduce_matrix(float** matrix, float** matrix_original, int nTypes)
{
	for(int i=0; i<nTypes; i++)
	{
		for(int j=0; j < nTypes; j++)
		{
			matrix_original[i][j] = matrix[i][j];
		}
	}

	for(int i=0; i<nTypes; i++)
	{
		float minimal_value = 1000000;
		for(int j=0; j < nTypes; j++)
		{
			if(matrix[i][j] < minimal_value)
			{
				minimal_value = matrix[i][j];
			}
		}
		for(int j =0; j < nTypes; j++)
		{
			matrix[i][j] -= minimal_value;
		}
	}

	for(int j=0; j<nTypes; j++)
	{
		float minimal_value = 1000000;
		for(int i =0; i<nTypes; i++)
		{
			if(matrix[i][j] < minimal_value)
			{
				minimal_value = matrix[i][j];
			}
		}
		for(int i=0; i<nTypes;i++)
		{
			matrix[i][j] -= minimal_value;
		}
	}
	return;
}