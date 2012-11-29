#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE buildlist, creates the list of atoms of residue rnum whose cart.
 * coordinates have to be rebuild after dihedral bnum is changed. Returns the 
 * vector with the lesser number of atoms that change orientation for ligands. 
 * In the case of amino acids the list contains the sidechain atoms that have 
 * to be rebuild. This vector contains the atom-number of those atoms sorted 
 * in a way that an atom always can be rebuild if the preceeding atoms were 
 * already build. In the case of ligand, it englofes everything covalently 
 * bonded to it.
 *****************************************************************************/
void buildlist(FA_Global* FA,atom* atoms,resid* residue,int rnum, int bnum, int *tot, int* lout){
	int a,b,i,j,k,l,flag;
	int num=0;
	int numa=0;
	int gpflag;
	int list[MAX_ATM_HET];
	int lista[MAX_ATM_HET];
	int rot;

	*tot=0;
	rot=residue[rnum].rot;
	//printf("rot: %d res.rot: %d\n", rot, residue[rnum].rot);
	j=residue[rnum].latm[rot]-residue[rnum].fatm[rot]+1;
	//printf("residue of rot: %d\tfatm: %d\tlatm: %d\n", rot, residue[rnum].fatm[0], residue[rnum].latm[0]);
	for(i=0;i<=j;i++){
		list[i]=0;
		lista[i]=0;
	}
	/* residue.type=0 -> amino acid 
	   residue.type=1 -> ligand     */
  
	if(residue[rnum].type==1){
		if(bnum != 0){
			a=atoms[residue[rnum].bond[bnum]].rec[0];
			b=atoms[residue[rnum].bond[bnum]].rec[1];
      
			//printf("a=%d %s, b=%d %s\n",a,atoms[a].name,b,atoms[b].name);
      
			/* adds to list the neighbours of a, exept b */
			for(i=1;i<=atoms[a].bond[0];i++){
				if(atoms[a].bond[i] != b){
					list[num]=atoms[a].bond[i];
					num++;
				}
			}
			/*num--;*/
			/* for each entry in list, takes their neighbours and if they 
			   are not a or b and are not already present in the list 
			   includes them into it
			*/
			for(i=0;i<=num;i++){
				for(j=1;j<=atoms[list[i]].bond[0];j++){
					k=atoms[list[i]].bond[j];
					if(k != a && k != b){
						flag=0;
						l=0;
						while(flag==0 && l<=num){
							if(k == list[l]){flag=1;}
							l++;
						}
						if(flag==0){
							list[num]=k;
							num++;
						}
					}
				}
			}

			/*for(k=0;k<num;k++){
			  printf("list[%d]=%d\n",k,list[k]);
			  }
			  PAUSE;
			*/
			/* check if the list contains the GPA atoms. if not it's ok otherwise the list
			   has to be inverted */
      
			gpflag=0;
			for(j=0;j<=2;j++){
				for(i=0;i<= num;i++){
					if(list[i]==residue[rnum].gpa[j]){
						gpflag=1;
						break;
					}
				}
			}

			/*printf("gpflag=%d\n",gpflag);
			  PAUSE;*/
      
			if(gpflag==1){
				rot=residue[rnum].rot;
				for(k=residue[rnum].fatm[rot];k<=residue[rnum].latm[rot];k++){
					flag=0;
					l=0;
					while(flag==0 && l<num){
						if(k == list[l]){flag=1;}
						l++;
					}
					if(flag==0){
						if(k!=a && k!=b){
							lista[numa]=k;
							numa++;
						}
					}
				}

				/*for(i=0;i<numa;i++){
				  printf("lista[%d]=%d\n",i,lista[i]);
				  }
				  PAUSE;*/

				for(i=0;i<numa;i++){list[i]=lista[i];}
				for(i=numa+1;i<=num;i++){list[i]=0;}
				num=numa;
			}

		}else{
			/* case the list is of all atoms for rigid body rotation & translation */
			num=0;
			//rot=residue[rnum].rot;
			//printf("rot: %d res.rot: %d\n", rot, residue[rnum].rot);
			//printf("Residue %d\tfatm: %d latm: %d\n", rnum, residue[rnum].fatm[0], residue[rnum].latm[0]);
			for(i=residue[rnum].fatm[rot];i<=residue[rnum].latm[rot];i++)
				list[num++]=i;

			//for (i=0;i<num;i++) {printf("list[%d]: %d\n",i,list[i]);}
			for(j=0;j<=2;j++){
				i=0;
				while(i<num){
					if(list[i] == residue[rnum].gpa[j]){break;}
					i++;
				}
				list[i]=list[j];
				list[j]=residue[rnum].gpa[j];
			}
		}

		//for(k=0;k<num;k++){printf("middle list[%d]=%d\n",k,list[k]);}
		//PAUSE;

		/* sorting of list according to reconstruction priority order */

		if(bnum == 0){                // the first three are ready to expert
			for(i=0;i<=2;i++){
				lout[*tot]=list[0];      // copy to final list the first
				(*tot)++;
				for(j=0;j<=num-1;j++){
					list[j]=list[j+1];   // copy down from list
				}
				num--;
			}
		}

		//for(k=0;k<num;k++){printf("after cuts list[%d]=%d\n",k,list[k]);}
		//printf("num=%d\n",num);
		//PAUSE;

		while(num > 0){
			//printf("here\n");
			for(i=0;i<num;i++){
				flag=0;
				//printf("list[%d]=%d\n",i,list[i]);
				for(j=0;j<=2;j++){
					l=0;
					while(flag==0 && l<num){
						//printf("atoms[%d].rec[%d]=%d list[%d]=%d\n",list[i],j,atoms[list[i]].rec[j],l,list[l]);
						//PAUSE;
						if(atoms[list[i]].rec[j] == list[l]){
							flag=1;
						}
						l++;
					}
					if(flag==1){break;}
				}

				if(flag==0){
					lout[*tot]=list[i];
					(*tot)++;
					for(j=i;j<num;j++){
						list[j]=list[j+1];
					}
					num--;
					//for(k=0;k<*tot;k++){printf("lout[%d]=%d\n",k,lout[k]);}
					//PAUSE;
				}
			}
		}

    
		for(k=0;k<*tot;k++){printf("lout[%d]= %d\n",k,atoms[lout[k]].number);}
		/*
		  for(k=0;k<*tot;k++){printf("lout[%d]=%d\n",k,lout[k]);}
		  PAUSE;
		*/
	}else{
		/* case of amino acids */

		/*printf("residue %d %s bond %d: %d\n",rnum,residue[rnum].name,
		  bnum,residue[rnum].bond[bnum]);
		  PAUSE;
		*/
		a=atoms[residue[rnum].bond[bnum]].rec[0];
		if(atoms[residue[rnum].bond[bnum]].rec[1] > a){
			a=atoms[residue[rnum].bond[bnum]].rec[0];
		}
		/*printf("a=%d\n",a);
		  PAUSE;*/
		rot=residue[rnum].rot;
		for(i=a+1;i<=residue[rnum].latm[rot];i++){
			/*printf("i=%d\n",i);
			  PAUSE;*/
			lout[num]=i;
			num++;
		}
		*tot=num;
	}
	return;
}
