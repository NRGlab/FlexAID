#include "flexaid.h"
#include "boinc.h"

void read_grid(FA_Global* FA,gridpoint** cleftgrid,char filename[]){

	FILE* file_ptr;  // read file stream
	char buffer[81]; // buffer to read line
	char field[7];   // 6 first character of the line
	char number[6];  // stores atom number
	char coor[9];    // stores coordinates
	int i,j,k;       // counters
	int flag;        // min flag
	float dis,mindis;    // minimal distance between intersections (spacer length)

	//printf("allocating memory for cleftgrid\n");
	(*cleftgrid) = (gridpoint*)malloc(FA->MIN_CLEFTGRID_POINTS*sizeof(gridpoint));    
  
	if ((*cleftgrid) == NULL){
		fprintf(stderr,"ERROR: memory allocation error for cleftgrid\n");
		Terminate(2);
	}
  
	//memset((*cleftgrid),0,FA->MIN_CLEFTGRID_POINTS*sizeof(gridpoint));

	file_ptr=NULL;
	if (!OpenFile_B(filename,"r",&file_ptr))
		Terminate(8);
  
	FA->num_grd = 1; // set counter to 1 because 0 is the INI conformation of the ligand
  
	while (fgets(buffer, sizeof(buffer), file_ptr)!=NULL){    
    
		if (FA->num_grd==FA->MIN_CLEFTGRID_POINTS){
			FA->MIN_CLEFTGRID_POINTS *= 2;
      
			//printf("re-allocating memory for cleftgrid\n");
			(*cleftgrid) = (gridpoint*)realloc((*cleftgrid),FA->MIN_CLEFTGRID_POINTS*sizeof(gridpoint));
			if ((*cleftgrid) == NULL){
				fprintf(stderr,"ERROR: memory reallocation error for cleftgrid\n");
				Terminate(2);
			}

			//memset(&(*cleftgrid)[FA->MIN_CLEFTGRID_POINTS/2],0,FA->MIN_CLEFTGRID_POINTS/2*sizeof(gridpoint));
		}
      

		for (i=0;i<6;i++) {
			field[i] = buffer[i];
		}
		field[6] = '\0';
    
		if (strcmp(field, "ATOM  ") == 0){
      
			for (i=6;i<11;i++) {
				number[i-6] = buffer[i];
			}
			number[5] = '\0';
			sscanf(number, "%d", &(*cleftgrid)[FA->num_grd].number);
      
			for (j=0;j<=2;j++){
				k=0;
				for (i=30+j*8;i<30+(j+1)*8;i++){
					coor[k] = buffer[i];
					k++;
				}
				coor[8] = '\0';
				sscanf(coor, "%f", &(*cleftgrid)[FA->num_grd].coor[j]);
			}

			FA->num_grd++;
      
		}
	}
  
	printf("num_grd: %d\n", FA->num_grd);

	CloseFile_B(&file_ptr,"r");

  
	flag=1;
	mindis=1000.0f;
	//Calculate spacer length
	for(i=2;i<FA->num_grd;i++){
		dis=distance((*cleftgrid)[1].coor,(*cleftgrid)[i].coor);
		if (flag) {
			mindis=dis;
			flag=0;
		}else{
			if (dis < mindis) mindis=dis;
		}
	}

	FA->spacer_length=mindis;
  
	return;
}
