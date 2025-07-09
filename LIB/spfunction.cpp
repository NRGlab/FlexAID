#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE cffunction calculates CF for a residue.
 * a residue may be either a protein amino acid or a hetero group.
 ******************************************************************************/
int spfunction(FA_Global* FA,atom* atoms,resid* residue){
  int    i,j,l,k,m,n,a,b,c;                    // dumb counters REMOVED s,p,k 

  int    nbonded;
  float  dis;

  cfstr* cfs;
  OptRes* optres;
  constraint* cons;

  float ang;
  int   npoints;
  double constant;

  float radA,radB,radC;
  float rAB;
  int   contnum;

  int    type;
  int    fatm;
  int    covalent;

  //int   blist[20];
  //int   nblist=0;
  // int   bloops;

  float  spoints[MAX_SPHERE_POINTS][3];
  float  spoints_flag[MAX_SPHERE_POINTS]; // flag indicating if sphere point was chosen

  float  Kwall;
  float  Kdist;
  float  Kangle;

  Kwall=1.0e6;
  Kdist=1.0e3;
  Kangle=1.0e2;
  
  for(i=0;i<FA->num_optres;i++){
    FA->optres[i].cf.rclash = 0;
    FA->optres[i].cf.com = 0.0;
    FA->optres[i].cf.con = 0.0;
    FA->optres[i].cf.sas = 0.0;
    FA->optres[i].cf.wal = 0.0;    
  }

  for(i=1;i<=FA->res_cnt;i++){
    
    optres = NULL;
    cfs = NULL;
    
    for(j=residue[i].fatm[residue[i].rot];j<=residue[i].latm[residue[i].rot];j++){
      
      if(atoms[j].optres != NULL){
	
	optres = atoms[j].optres;
	type = optres->type;
	cfs = &optres->cf;
	
      }else{

	// non-optimizable residue
	continue;
      }
      
      if(type == 0 && atoms[j].isbb){continue;}

      for(c=0;c<MAX_SPHERE_POINTS;c++){
	spoints_flag[c]=0;
      }

      radA = atoms[j].radius;

      // approximation of area of one point of sphere
      constant=4.0*PI*pow((double)(radA+Rw),2.0)/(double)FA->tspoints;
      
      contnum = 0;
      nbonded = 0;

      
      if(atoms[j].ncons > 0){
	for(m=0;m<atoms[j].ncons;m++){

	  radC = atoms[j].number==atoms[j].cons[m]->anum1?
	    atoms[FA->num_atm[atoms[j].cons[m]->anum2]].radius:
	    atoms[FA->num_atm[atoms[j].cons[m]->anum1]].radius;

	  // maximum penalty value (starting penalty)
	  // default value if atoms are not interacting
	  cfs->con += Kangle*(radA+radC+2.0f*Rw);
	  cfs->con += Kdist*(radA+radC+2.0f*Rw);
	  
	}
      }

      // find contacts
      for(k=1;k<=FA->res_cnt;k++){
	
	for(l=residue[k].fatm[residue[k].rot];l<=residue[k].latm[residue[k].rot];l++){
	  
	  // intramolecular atoms
	  if(k==i){
	    fatm = residue[i].fatm[residue[i].rot];
	    if(residue[i].bonded[j-fatm][l-fatm] >= 0){
	      // atom is bonded
	      nbonded++;
	      continue;
	    }
	  }
	  	  
	  radB=atoms[l].radius;
	  rAB=radA+radB;
	  
	  dis=sqrdist(atoms[j].coor,atoms[l].coor);
	  
	  // is atom in contact range.
	  // if NO, discard atom l
	  if(dis > (rAB+2.0*Rw)*(rAB+2.0*Rw)){continue;}
	  
	  covalent = 0;
	  cons=NULL;

	  if(atoms[j].ncons > 0 && atoms[l].ncons > 0){
	    for(m=0;m<atoms[j].ncons;m++){
	      for(n=0;n<atoms[l].ncons;n++){
		if(atoms[j].cons[m]->id == atoms[l].cons[n]->id){
		  cons = atoms[j].cons[m];
		  break;
		}
	      }
	      if(cons != NULL){break;}
	    }	 
	    
	    if(cons != NULL){
	      
	      // covalent constraint
	      if(cons->type == 1){
		covalent = 1;
				
		ang = angle(atoms[j].coor,atoms[l].coor,atoms[atoms[l].bond[1]].coor);
		cfs->con -= Kangle*(rAB+2.0f*Rw)*(1.0f-(fabs(ang-120.0f)/120.0f));

		cfs->con -= Kdist*(rAB+2.0f*Rw)*(1.0f-(fabs(dis-1.5f)/1.5f));
	
	      }

	      // interaction constraint
	      else{
		

	      }
	      
	    } // end of constraint found
	    
	  } // end of constraint
       	
	  // WALL term
	  if(dis <= (FA->permeability*rAB)*(FA->permeability*rAB)  &&
	     !covalent || dis <= 1.5*1.5){
	    
	    cfs->wal += Kwall*(pow(dis,-6.0f)-pow(FA->permeability*rAB,-12.0f));

	    if(dis <= (FA->dee_clash*rAB)){cfs->rclash=1;}

	  }

	  npoints=0;
	  // calculate complementarity using spheres of points
	  for(a=0;a<FA->tspoints;a++){
	    for(b=0;b<3;b++){
	      spoints[a][b] = FA->sphere[a][b]*(radA+Rw)+atoms[j].coor[b];
	    }

	    if(sqrdist(spoints[a],atoms[l].coor) <= (atoms[l].radius+Rw)*(atoms[l].radius+Rw)){
	      spoints_flag[a]=1;
	      npoints++;
	    }

	  }

	  cfs->com += FA->energy[atoms[j].type][atoms[l].type]*constant*npoints;


	  // number of contacts for atom j
	  contnum++;

	} // end of atom l
	
      } // end of residue k
      
      npoints=0;

      for(a=0;a<FA->tspoints;a++){
	if(spoints_flag[a]){npoints++;}
      }
      
      if(nbonded){ npoints += FA->tspoints/4*nbonded; }

      if(npoints > FA->tspoints){ npoints = FA->tspoints; }

      printf("number of points in contact for atom[%d] = %d\n",j,npoints);

      cfs->sas += (FA->tspoints-npoints)*constant;

    } // end of atom j
    
  } // end of residue i
  
  return(0);

} 
