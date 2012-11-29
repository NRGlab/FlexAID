#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE cffunction calculates CF for a residue.
 * a residue may be either a protein amino acid or a hetero group.
 ******************************************************************************/
cfstr cffunction(FA_Global* FA,atom* atoms,resid* residue,int num){
  int   i,j,l,k,m,n,b,c;                    // dumb counters REMOVED s,p,k 
  int   bond_flag;                          // flag to test bonded atoms 
  int   cont_flag;                          // flag to test atoms in contact 
  float dis;
  float distances[100];
  int   close2atom[100];
  int   num_close;
  int   rot;
  
  //int   blist[20];
  //int   nblist=0;
  // int   bloops;
  float spoint[MAX_SPHERE_POINTS][3];
  float constant;
  float rij;
  float Kwall=0.0;
  float Ewall=0.0;
  char  rangeswitch;
  cfstr cf;
  float com=0.0; 
  float nor=0.0;
  float overlap;

  //comment the following if not for full list of contacts.
  float freesurf=0.0;
  //float occ[100];
  //float com_res=0.0; 
  //int ncom;
  //end comment

  Kwall=1.0e6;
  overlap=0.9;
  cf.rclash=0;

  // type of calculation, a=includes intramolecular interactions
  //                      b=excludes intramolecular interactions  
  // a for flexible - b is equal to fortran version
  rangeswitch='a';

  // number of bonds away an atom is added to blist for the purpose
  // of checking the wall term
  // bloops=2;

  for(i=residue[num].fatm[residue[num].rot];i<=residue[num].latm[residue[num].rot];i++){

    constant=4.0*PI*pow((atoms[i].radius+Rw),2.0)/(float)FA->tspoints;
    
    // bondedlist(i,bloops,&nblist,blist);
    for (j=0;j<FA->num_het_atm;j++){
      if(FA->bondlist[j].num==i) {
	b=j;
	break;
      }
    }  

    // creates list of atoms j in contact with i sorted in increasing 
    // distance from i. keeps up to 100 closest atoms in contact.
    num_close=0;
    for(l=0;l<100;l++){close2atom[l]=0;};
    for(k=1;k<=FA->res_cnt;k++){
      if (residue[k].type==1 && k != num) { continue; }
      rot=residue[k].rot;
      for(j=residue[k].fatm[rot];j<=residue[k].latm[rot];j++){
	if(j != i){
	  dis=sqrdist(atoms[i].coor,atoms[j].coor);
	  rij=atoms[i].radius+atoms[j].radius;
	  if(dis <= (rij+2.0*Rw)*(rij+2.0*Rw)){ 
	    if(num_close == 0){
	      // No need to sort because 0 element in array
	      distances[num_close]=dis;
	      close2atom[num_close]=j;
	    }else{
	      // Sort if more than 0 elements in array
	      m=num_close-1;
	      while(distances[m]>dis && m >= 0){                
		distances[m+1]=distances[m];
		close2atom[m+1]=close2atom[m];
		m--;
	      }
	      distances[m+1]=dis;
	      close2atom[m+1]=j;
	    }
	    num_close++;
	    
	    // Considered as bonded only if its bloops bonds away from atom j
	    bond_flag=0;
	    for (c=0;c<FA->bondlist[b].tot;c++){
	      if (j==FA->bondlist[b].nbr[c]){
		bond_flag=1;
		break;
	      }
	    }
	    /*
	    for(l=1;l<nblist;l++){
	      if(blist[l] == j){
		bond_flag = 1;
		break;
	      }
	    }	 
	    */

	    if(bond_flag==0){
	      // calculates the wall term for overlaps
	      // If (Rij < R0) where R0=0.9(Ri+Rj)
	      if(dis <= (overlap*rij)*(overlap*rij)){
		//printf("Atom overlap between [%d]-[%d]\n",i,j);
		Ewall += Kwall*(pow(dis,-6.0)-pow(0.9*rij,-12.0));
	      }	      
	      if(dis <= 4.0){
		// Intramolecular clash because of flexibility
		if (k==num){cf.rclash=1;}
		// Extramolecular clash with a rigid side-chain from binding pocket
		else{if (atoms[j].recs=='r'){cf.rclash=1;}}
	      }
	    }
	  }
	}
      }
    }
    
    n=0;
    for(l=0;l<FA->tspoints;l++){
      //create sphere point l around atom i
      for(m=0;m<=2;m++){
	spoint[l][m]=(atoms[i].radius+Rw)*FA->sphere[l][m]+atoms[i].coor[m];
      }

      //check if point l is not in contact to residue num
      cont_flag=0;
      rot=residue[num].rot;
      for(j=residue[num].fatm[rot];j<=residue[num].latm[rot];j++){
	if(j != i){
	  if(sqrdist(spoint[l],atoms[j].coor) <= (atoms[j].radius+Rw)*(atoms[j].radius+Rw)){
	    cont_flag=1;
	    break;
	  }
	}
      }
      if(cont_flag==0){n++;}
    }
    nor += (float)n*constant;


    //comment the following if not for full list of contacts.
    //for(m=0;m<=num_close;m++){occ[m]=0.0;}
    //end comment

    //printf("num_close: %d\n",num_close);
    for (l=0;l<FA->tspoints;l++){
      m=0;
      while(m<num_close){
	dis=sqrdist(spoint[l],atoms[close2atom[m]].coor);
	if(dis <= (atoms[close2atom[m]].radius + Rw)*(atoms[close2atom[m]].radius + Rw)){
	  bond_flag=0;
	  
	  for(n=1;n<FA->bondlist[c].tot;n++){
	    if(FA->bondlist[c].nbr[n] == close2atom[m]){
	      bond_flag=1;
	      break;
	    }
	  }

	  // excludes intramolecular interactions
	  if(rangeswitch == 'b' && bond_flag == 0){
	    if(atoms[close2atom[m]].ofres == num){
	      bond_flag=1;
	    }
	  }
	  break;
	}
	m++;
      }
      if(m != num_close){
	if(bond_flag == 0){
	  //comment the following if not for full list of contacts.
	  //occ[m] += constant;
	  //end comment
	  com += FA->energy[atoms[i].type][atoms[close2atom[m]].type]*constant;
	}
      }
      //comment the following if not for full list of contacts.
      else{freesurf += constant;}
      //end comment
    }
    
    //comment the following if not for full list of contacts.
    /*
      printf("Atom %4d %s (type %d) has contacts with:\n",atoms[i].number, 
      atoms[i].name,atoms[i].type);
      for(m=0;m<=num_close;m++){
      if(occ[m] != 0.0){
      com_res += FA->energy[atoms[i].type][atoms[close2atom[m]].type]*occ[m];
      printf("%s %3d %s (type %d->%4.2f) SC=%6.3f COM=%8.3f\n",
      residue[atoms[close2atom[m]].ofres].name,
      residue[atoms[close2atom[m]].ofres].number,
      atoms[close2atom[m]].name,atoms[close2atom[m]].type,
      FA->energy[atoms[i].type][atoms[close2atom[m]].type],occ[m],com_res);
      }
      }
      PAUSE;
    */
    //end comment

  }
  //comment the following if not for full list of contacts.
  //printf("COM=%f NOR=%f Ewall=%f free=%f\n",com,nor,Ewall,freesurf);
  //end comment
  //com = (com - Ewall)/ncom ;

  cf.com=com-2.0*freesurf;
  //cf.com=com;
  cf.nor=nor;
  cf.wal=Ewall;

  return(cf);
} 
