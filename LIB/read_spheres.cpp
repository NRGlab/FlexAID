#include "flexaid.h"
#include "boinc.h"

sphere* read_spheres(char filename[]){

	FILE* file_ptr;             // read file stream
	char buffer[81];            // buffer to read line
	char field[7];              // 6 first character of the line
	char coor[9];               // coordinate
	char radius[6];             // radius
	int i,j,k;                  // dummy counters
	sphere* spheres = NULL;     // spheres list
	
	file_ptr=NULL;
	if (!OpenFile_B(filename,"r",&file_ptr))
		Terminate(8);
  
	while (fgets(buffer, sizeof(buffer), file_ptr) != NULL){
		
		//0         1         2         3         4         5         6         7
		//0123456789012345678901234567890123456789012345678901234567890123456789
		//ATOM   1102  C   SPH Z   1      34.069  28.877   7.194  1.00  1.51 
		for (i=0;i<6;i++) field[i] = buffer[i];
		field[6] = '\0';
    
		if (strcmp(field, "ATOM  ") == 0){
			
			sphere* _sphere;
			_sphere = (sphere*)malloc(sizeof(sphere));
			if(_sphere == NULL){
				fprintf(stderr,"ERROR: memory allocation error for spheres (LOCCLF)\n");
				Terminate(2);
			}
			
			
			for (j=0;j<=2;j++){
				k=0;
				for (i=30+j*8;i<30+(j+1)*8;i++)
					coor[k++] = buffer[i];
				coor[8] = '\0';
				
				sscanf(coor, "%f", &_sphere->center[j]);
			}
			
			for(i=0; i<5; i++)
				radius[i] = buffer[i+61];
			radius[5] = '\0';

			sscanf(radius, "%f", &_sphere->radius);
			
			/*
			printf("new sphere\nradius: %f\ncoor[0]: %.3f coor[1]: %.3f coor[2]: %.3f\n",
			       _sphere->radius, _sphere->center[0], _sphere->center[1], _sphere->center[2]);
			*/

			_sphere->prev = spheres;
			spheres = _sphere;
		}
	}
  
	return spheres;
}
