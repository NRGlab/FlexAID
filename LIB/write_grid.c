#include "gaboom.h"
#include "boinc.h"

/************************************************************************/
/*******   THIS FUNCTION IS USED TO WRITE TO A FILE THE LINES ***********/
/*******   THAT ARE USED TO VISUALIZE THE GRID IN pyMol       ***********/
/*******   PYTHON SCRIPT grid.py IS NECESSARY IN ORDER TO VIEW***********/
/*******   FORMAT OF A LINE: X1 Y1 Z1 X2 Y2 Z2\n              ***********/
/************************************************************************/

void write_grid(FA_Global* FA,const gridpoint* cleftgrid,char gridfilename[]) {


	FILE *outfile_ptr;

	const char chn[10] = { 'A', 'B', 'C', 'D', 'E',
			       'F', 'G', 'H', 'I', 'J' };
	
	outfile_ptr=NULL;
	if (!OpenFile_B(gridfilename,"w",&outfile_ptr)){
		Terminate(6);
	}else{
		//printf("spacer_length=%.3f\n",FA->spacer_length);
		
		int i,j;
		int c = 0;

		for(i=1,j=1;i<FA->num_grd;i++,j++){

			if(i % 100000 == 0){
				c++;
				j=1;
			}

			fprintf(outfile_ptr, 
				"ATOM  %5d  C   GRD %c   1    %8.3f%8.3f%8.3f  0.00  0.00\n", 
				j, chn[c], cleftgrid[i].coor[0], cleftgrid[i].coor[1], cleftgrid[i].coor[2]);
			
		}
		CloseFile_B(&outfile_ptr,"w");
	}

	return;
}
