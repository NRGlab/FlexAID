#include "Vcontacts.h"

#define DEBUG_LEVEL 0

int vcfunction(FA_Global* FA,VC_Global* VC,atom* atoms,resid* residue)
{
	int    i,j,k;
	double SAS,area;
	cfstr* cfs = NULL;
	int    currindex;
	int    contnum;       // number of contacts (excluding bloops away atoms)
	//float *coorA, *coorB; // coordinate pointers
	double  radA,radoA, radB, radC, rAB;
	double  permea,dee_clash;

	//float  Ewall;
	//float  Ewall_atm=0;

	//float  com_atm;
	//int    rot=0;
	int    rnum=0;
	int    type=1;
	int    fatm=0;

	double dist_opt = 0.0;
	double complementarity;
	//int   intra;
	//int   intraf;

	int atomzero=0;    // atom index in atoms structure of source atom
	int atomcont=0;    // atom index in ****************** contacting atom

	int nconszero=0; 
	int nconscont=0;
	int covalent=0;
	constraint* cons = NULL;

#if DEBUG_LEVEL > 0
	cfstr cf_atom;
#endif
  
	//float  matrix[NTYPES+1][NTYPES+1];

	// reset all values pointed
	memset(VC->Calc,0,FA->atm_cnt_real*sizeof(atomsas));
	memset(VC->Calclist,0,FA->atm_cnt_real*sizeof(int));
	memset(VC->ca_index,0,FA->atm_cnt_real*sizeof(int));
	memset(VC->seed,0,3*FA->atm_cnt_real*sizeof(int));
	memset(VC->contlist,0,10000*sizeof(contactlist));
	memset(FA->contacts,0,100000*sizeof(int));

	// reset CF values
	for(j=0; j<FA->num_optres; ++j){
		FA->optres[j].cf.rclash=0;
		FA->optres[j].cf.wal=0.0;
		FA->optres[j].cf.com=0.0;
		FA->optres[j].cf.con=0.0;
		FA->optres[j].cf.sas=0.0;
	}
    
	permea = (double)FA->permeability;
	dee_clash = (double)FA->dee_clash;
	
	/*
	// empty surface matrix values
	for(i=0;i<=NTYPES;++i){
	for(j=0;j<=NTYPES;++j){
	matrix[i][j]=0.0;
	}
	}
	*/  
    
	//printf("=============NEW INDIVIDUAL==============\n");
    
    
	if(Vcontacts(FA,atoms,residue,VC) == -1){

        FA->skipped++;

		free(VC->ca_rec);
		free(VC->box);

		return(1);
	}
    
	for(i=0; i<FA->atm_cnt_real; ++i) {
    
		cfs = NULL;

		// atom from which contacts are calculated
		atomzero = FA->num_atm[VC->Calc[i].atomnum];

		// number of constraints for atomzero
		nconszero=atoms[atomzero].ncons;

		if(atoms[atomzero].optres != NULL)
		{
			// the residue optimizable
			rnum = atoms[atomzero].optres->rnum;
			type = atoms[atomzero].optres->type;
			cfs = &atoms[atomzero].optres->cf;
		}
		else
		{
			continue;
		}

		if(type == 0 && atoms[atomzero].isbb)
		{
			continue;
		}

		//printf("-------------------------------\nAtom[%4d]=%s in Residue[%4d]\n",VC->Calc[i].atomnum,atoms[FA->num_atm[VC->Calc[i].atomnum]].name,VC->Calc[i].resnum);
    
    
		//com_atm=0.0;
		//Ewall_atm=0.0;

		radA  = (double)VC->Calc[i].radius;
		radoA = radA + Rw;

		SAS = 4.0*PI*radoA*radoA;

#if DEBUG_LEVEL > 0
		cf_atom.sas = 0.0;
#endif

		if(atoms[atomzero].ncons > 0){

			for(j=0;j<atoms[atomzero].ncons;j++){
	
				radC = atoms[atomzero].number==atoms[atomzero].cons[j]->anum1 ? 
					(double)atoms[FA->num_atm[atoms[atomzero].cons[j]->anum2]].radius:
					(double)atoms[FA->num_atm[atoms[atomzero].cons[j]->anum1]].radius;
	
				// maximum penalty value (starting penalty)
				// default value if atoms are not interacting
				//cfs->con += KANGLE*(radA+radC+2.0*Rw);

				if(atoms[atomzero].cons[j]->type == 1){
					cfs->con += KDIST;
				}
				
				//printf("constraint for atom[%d]: %.3f\n", atoms[atomzero].number,cfs->con);
			}
      
			/*
			  printf("default constraint value: %.3f\n",cfs->con);
			  getchar();
			*/
		}


#if DEBUG_LEVEL > 0
		printf("==================================================================================\n");
		printf("ATOM :: RES C RNUM  ANUM  T  RAD ::   COMPL   DIST   AREA ::     CF.COM     CF.WAL\n");
		printf("----------------------------------------------------------------------------------\n");
		printf("ATOM :: %3s %c %4d %5d %2d %4.2f\n",
		       residue[atoms[atomzero].ofres].name,
		       residue[atoms[atomzero].ofres].chn,
		       residue[atoms[atomzero].ofres].number,
		       atoms[atomzero].number,
		       atoms[atomzero].type,
		       atoms[atomzero].radius);
#endif


		currindex = VC->ca_index[i];

		contnum = 0;

		while(currindex != -1) {

			//printf("[%6.3f] Contact between atom[%5d]-T=[%d]---atom[%5d]T=[%d] with AREA[%8.3f] DIST[%4.3f]\n",FA->energy[VC->Calc[i].type][VC->Calc[VC->ca_rec[currindex].atom].type],VC->Calc[i].atomnum,VC->Calc[i].type,VC->Calc[VC->ca_rec[currindex].atom].atomnum,VC->Calc[VC->ca_rec[currindex].atom].type,VC->ca_rec[currindex].area,VC->ca_rec[currindex].dist); 

			complementarity = 0.0;
			area = VC->ca_rec[currindex].area;

			SAS -= area;

			// number of contacts counter
			contnum++;

#if DEBUG_LEVEL > 0
			cf_atom.com  =  0.0;
			cf_atom.wal  =  0.0;
#endif
			
			/*
			if ( VC->ca_rec[currindex].dist < 0.0001 ) { 
				currindex = VC->ca_rec[currindex].prev;
				continue; 
			}
			*/

			// atom in contact with atom zero
			atomcont = FA->num_atm[VC->Calc[VC->ca_rec[currindex].atom].atomnum];

			// is contact atom bonded to atom zero
			// if YES, skip contact atom
			if(atoms[atomcont].ofres == atoms[atomzero].ofres)
			{
	  
				// get first atom of residue
				fatm = residue[rnum].fatm[0];

				if(residue[rnum].bonded[atomcont-fatm][atomzero-fatm] >= 0)
				{

#if DEBUG_LEVEL > 2
					printf("    (B) %3s %c %4d %5d %2d %4.2f :: %7.4f %6.2f %6.2f :: %10.3f %10.3f\n",
					       VC->Calc[VC->ca_rec[currindex].atom].res,
					       VC->Calc[VC->ca_rec[currindex].atom].chn,
					       VC->Calc[VC->ca_rec[currindex].atom].resnum,
					       VC->Calc[VC->ca_rec[currindex].atom].atomnum,
					       VC->Calc[VC->ca_rec[currindex].atom].type,
					       VC->Calc[VC->ca_rec[currindex].atom].radius,
					       
					       complementarity,
					       VC->ca_rec[currindex].dist,
					       VC->ca_rec[currindex].area,
					       
					       cf_atom.com, cf_atom.wal);

#endif
					currindex = VC->ca_rec[currindex].prev;
					continue;	  
				}
			}

			if(FA->contacts[VC->Calc[VC->ca_rec[currindex].atom].atomnum]){
				//printf("%d already calculated\n",VC->Calc[VC->ca_rec[currindex].atom].atomnum );
				currindex = VC->ca_rec[currindex].prev;
				continue;
			}

			// covalently bonded flag
			covalent = 0;

			// number of constraints for contact atom
			nconscont=atoms[atomcont].ncons;
      
			// complementarity value between the 2 atom types
			if(VC->Calc[i].type == 0 || VC->Calc[VC->ca_rec[currindex].atom].type == 0){
				complementarity = 0.0;
			}else{
				complementarity = (double)FA->energy[VC->Calc[i].type][VC->Calc[VC->ca_rec[currindex].atom].type];
			}
      
			radB  = (double)VC->Calc[VC->ca_rec[currindex].atom].radius;
			rAB   = radA+radB;

			// do contacting atoms have the same constraint
			cons=NULL;
			if(nconszero > 0 && nconscont > 0){
				for(j=0;j<nconszero;j++){
					for(k=0;k<nconscont;k++){
						if(atoms[atomzero].cons[j]->id == atoms[atomcont].cons[k]->id){
							cons = &FA->constraints[atoms[atomzero].cons[j]->id];
							break;
						}
					}
					if(cons != NULL){break;}
				}
				
				if(cons != NULL){
					
					if(cons->type == 1){
						
						covalent=1;
						dist_opt = cons->bond_len;
						
						//ang = angle(atoms[atomzero].coor,atoms[atomcont].coor,atoms[atoms[atomcont].bond[1]].coor);
						//cfs->con -= KANGLE*cons->max_ang*GetValueFromGaussian(ang,120.0,cons->max_ang);
						
						//float gaus=GetValueFromGaussian(VC->ca_rec[currindex].dist,dist_opt,cons->max_dist);
						//printf("gaus: %.2f\n",gaus);
						//getchar();
						
						printf("cons->max_dist=%.3f - VC->ca_rec[currindex].dist=%.3f - dist_opt=%.3f - GetValueFromGaussian=%.3f\n", cons->max_dist, VC->ca_rec[currindex].dist, dist_opt,GetValueFromGaussian(VC->ca_rec[currindex].dist,dist_opt,cons->max_dist));
						
						printf("before: %.3f\n",cfs->con);
						cfs->con -= KDIST * GetValueFromGaussian(VC->ca_rec[currindex].dist,dist_opt,cons->max_dist);
						printf("after: %.3f\n",cfs->con);
						//getchar();
						
						/*
						  printf("removing from con: %.3f\n", //KANGLE*cons->max_ang*GetValueFromGaussian(ang,120.0,cons->max_ang));
						  KDIST*cons->max_dist*GetValueFromGaussian(VC->ca_rec[currindex].dist,dist_opt,cons->max_dist));
						  getchar();
						*/
						
					}else{

						// interaction constraint
						if(cons->force_interaction){
							complementarity = 
								complementarity < 0 ? 
								-1.0 * complementarity:
								complementarity;
						}

						complementarity *= cons->interaction_factor;
		  
						/*
						  printf("found interaction constraint: %.2f (%d)\n",
						  cons->interaction_factor,cons->force_interaction);
						*/
					}

				}
	
				//printf("constraint[%d] applies.\n",cons->id);
			}
            
			//	coorB = VC->Calc[VC->ca_rec[currindex].atom].coor;

			// CHECK IF CLASH
			float clash_distance = permea*rAB;
			if(covalent && permea*dist_opt < permea*rAB){
				clash_distance = permea*dist_opt;
			}
                        
			if (VC->ca_rec[currindex].dist < clash_distance){
				
				/*
				  if(covalent){
				  int atmi=FA->num_atm[VC->Calc[i].atomnum];
				  int atmj=FA->num_atm[VC->Calc[VC->ca_rec[currindex].atom].atomnum];
				  int resi=atoms[atmi].ofres;
				  int resj=atoms[atmj].ofres;
				  
				  printf("((Atom overlap between %s%d%c[%d](%s)-%s%d%c[%d](%s))) with DIST=%.3f with CLASH_DIST=%.3f\n",
				  residue[resi].name,residue[resi].number,residue[resi].chn,
				  atoms[atmi].number,atoms[atmi].name,
				  residue[resj].name,residue[resj].number,residue[resj].chn,
				  atoms[atmj].number,atoms[atmj].name,
				  VC->ca_rec[currindex].dist,clash_distance);
				  }
				*/
				
				//Ewall_atm += KWALL*(pow(VC->ca_rec[currindex].dist*VC->ca_rec[currindex].dist,-6.0)-pow(0.9*rAB,-12.0));
				cfs->wal += KWALL*(pow(VC->ca_rec[currindex].dist,-12.0)-pow(permea*rAB,-12.0));

#if DEBUG_LEVEL > 0
				cf_atom.wal += KWALL*(pow(VC->ca_rec[currindex].dist,-12.0)-pow(permea*rAB,-12.0));
#endif

				// Treat everything as rigid
				if ( VC->ca_rec[currindex].dist <= dee_clash*rAB ) { cfs->rclash=1; }

			}
            
			// ATOM COMPLEMENTARITY
			//com_atm += (double)FA->energy[VC->Calc[i].type][VC->Calc[VC->ca_rec[currindex].atom].type]*area;
			//cfs->com += intraf*complementarity*area;
			if( covalent == 0 ) { 
				cfs->com += complementarity * area;
			
#if DEBUG_LEVEL > 0	
				cf_atom.com += complementarity * area;
#endif

			}
      
			/*
			// generate surface area matrix
			matrix[VC->Calc[i].type][VC->Calc[VC->ca_rec[currindex].atom].type] += area;
      
			if(VC->Calc[i].type != VC->Calc[VC->ca_rec[currindex].atom].type){
			matrix[VC->Calc[VC->ca_rec[currindex].atom].type][VC->Calc[i].type] += area;
			}
			*/
      			
#if DEBUG_LEVEL > 0
			printf("        %3s %c %4d %5d %2d %4.2f :: %7.4f %6.2f %6.2f :: %10.3f %10.3f\n",
			       VC->Calc[VC->ca_rec[currindex].atom].res,
			       VC->Calc[VC->ca_rec[currindex].atom].chn,
			       VC->Calc[VC->ca_rec[currindex].atom].resnum,
			       VC->Calc[VC->ca_rec[currindex].atom].atomnum,
			       VC->Calc[VC->ca_rec[currindex].atom].type,
			       VC->Calc[VC->ca_rec[currindex].atom].radius,
			       
			       complementarity,
			       VC->ca_rec[currindex].dist,
			       VC->ca_rec[currindex].area,
			       
			       cf_atom.com, cf_atom.wal);
			
#endif

			// skip to next contact
			currindex = VC->ca_rec[currindex].prev;
		}    
        
		//    printf("Atom[%d]=%d has %d contacts\n",VC->Calc[i].atomnum,VC->Calc[i].atomnum,contnum);
		//    printf("Atom[%d] COM=[%8.2f]\tWAL=[%8.2f]\n",VC->Calc[i].atomnum,com_atm,Ewall_atm);
    
		if(SAS < 0.0){ SAS = 0.0; }

		if(FA->by_solventtype){
			cfs->sas += (double)FA->energy[VC->Calc[i].type][FA->by_solventtype] * SAS;
			
#if DEBUG_LEVEL > 0
			cf_atom.sas += (double)FA->energy[VC->Calc[i].type][FA->by_solventtype] * SAS;
#endif
		}else{
			cfs->sas += (double)FA->solventterm * SAS;

#if DEBUG_LEVEL > 0
			cf_atom.sas += (double)FA->solventterm * SAS;
#endif
		}

		
		FA->contacts[VC->Calc[i].atomnum] = 1;
		//printf("%d calculated\n", VC->Calc[i].atomnum);


#if DEBUG_LEVEL > 1
		printf("CF.SAS is %.3f (Area=%.3f) for %d contacts\n", 
		       cf_atom.sas, SAS, contnum); 
#endif

	}


	// Penalize Freesurf.
  
	//printf("FREESURF(SAS)=[%8.2f]\n",SAStot);
	//printf("FINAL SUM COM=[%8.2f]\tWAL=[%8.2f]\n",com,Ewall);
  
	//print_surfmat(matrix,"surf.mat");
#if DEBUG_LEVEL > 0
	printf("\n");
	printf("CF.sum = %.3f\n", cfs->com + cfs->sas + cfs->wal);
	printf("CF.com = %.3f\n", cfs->com);
	printf("CF.sas = %.3f\n", cfs->sas);
	printf("CF.wal = %.3f\n", cfs->wal);
	getchar(); 
#endif
  	
	free(VC->ca_rec);
	//printf("free-ing %p\n",VC->ca_rec);
  
	free(VC->box);

  
	return(0);
  
}
