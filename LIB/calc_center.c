#include "flexaid.h"
/***********************************************************
 * This subroutine calculates the coordinates of the center 
 * of geometry of all protein atoms.
 * It also finds the minimum and maximum coordinates
 * for X, Y and Z of all proteins atoms (including HET)
 ***********************************************************/
void calc_center(FA_Global* FA,atom* atoms,resid* residue){
	int i,j,k;
	int tot;
	//int atm;
	float d;
	int rot;

	for (j=0;j<3;j++){
		FA->globalmin[j]=atoms[1].coor[j];
		FA->globalmax[j]=atoms[1].coor[j];
	}

	for(j=0;j<=2;j++){
		tot=0;
		FA->ori[j] = 0.0f;
		for(k=1;k<=FA->res_cnt;k++){
			// Skip HET group atoms
			// if (residue[k].type==1) continue;
      
			//rot=residue[k].rot;
			for(i=residue[k].fatm[0];i<=residue[k].latm[0];i++){
				FA->ori[j] += atoms[i].coor[j];
				if (atoms[i].coor[j] < FA->globalmin[j]) FA->globalmin[j]=atoms[i].coor[j];
				if (atoms[i].coor[j] > FA->globalmax[j]) FA->globalmax[j]=atoms[i].coor[j];
				//printf("%7.3f\n",atoms[i].coor[j]);
				tot++;
			}
		}
		FA->ori[j] /= (float)tot;
	}

	FA->maxdst=0;
	for(k=1;k<=FA->res_cnt;k++){
		rot=residue[k].rot;
		for(i=residue[k].fatm[rot];i<=residue[k].latm[rot];i++){
			d=distance(FA->ori,atoms[i].coor);
			if(d >= FA->maxdst) FA->maxdst=d;
		}
	}

	FA->maxwidth=0.0f;
	for(j=0; j<3; j++) {
		if((FA->globalmax[j]-FA->globalmin[j]) > FA->maxwidth) {
			FA->maxwidth = (FA->globalmax[j]-FA->globalmin[j]);
		}
	}
	//FA->maxdst=(int)FA->maxdst+1.0;
  
	for(j=0;j<3;j++){
		printf("globalmin[%d]=%8.3f\tglobalmax[%d]=%8.3f\n",j,FA->globalmin[j],j,FA->globalmax[j]);
	}
	//PAUSE;

	return;
}
