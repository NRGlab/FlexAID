#include "Vcontacts.h"

/*#define DEBUG_LEVEL 0*/

double vcfunction(FA_Global* FA,VC_Global* VC,atom* atoms,resid* residue, std::vector< std::pair<int,int> > & intraclashes, bool* error)
{
	int    rnum=0;
	int    type=1;
	
	// reset all values pointed
	memset(FA->contacts,0,100000*sizeof(int));
	memset(FA->contributions,0.0f,FA->ntypes*FA->ntypes*sizeof(float));
	
	// reset CF values
	for(int j=0; j<FA->num_optres; ++j){
		FA->optres[j].cf.rclash=0;
		FA->optres[j].cf.wal=0.0;
		FA->optres[j].cf.com=0.0;
		FA->optres[j].cf.con=0.0;
		FA->optres[j].cf.sas=0.0;
		FA->optres[j].cf.ligsolv=0.0;
		FA->optres[j].cf.tarsolv=0.0;
		FA->optres[j].cf.totsas=0.0;
	}
	
	double permea = (double)FA->permeability;
	double dee_clash = (double)FA->dee_clash;
	
	// allocate
	//float  matrix[FA->ntypes*FA->ntypes];
	/*
	// empty surface matrix values
	for(i=0;i<FA->ntypes;++i){
	for(j=0;j<FA->ntypes;++j){
	matrix[i][j]=0.0;
	}
	}
	*/  
    
	//printf("=============NEW INDIVIDUAL==============\n");
    
	double clash_value;
	*error = false;
	int rv = Vcontacts(FA,atoms,residue,VC,&clash_value,false);
	if(rv < 0){
		*error = true;
		
		for(int i=0;i<FA->atm_cnt_real;i++){
			if(VC->Calc[i].score){
				VC->Calc[i].atom = NULL;
			}
		}
		
		if(!FA->vindex){ free(VC->box); }
		
		if(rv == -1){
			FA->skipped++;
			return(POLYHEDRON_PENALTY);
		}else if(rv == -2){
			FA->clashed++;
			return(clash_value);
		}
	}
	
	for(int i=0; i<FA->atm_cnt_real; ++i) {
		
		cfstr* cfs = NULL;
#if DEBUG_LEVEL > 0
		cfstr cfs_atom;
#endif

		// atom from which contacts are calculated
		int atomzero = FA->num_atm[VC->Calc[i].atom->number];
		
		// number of constraints for atomzero
		int nconszero = atoms[atomzero].ncons;
		
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
		
		//printf("-------------------------------\nAtom[%4d]=%s in Residue[%4d]\n",VC->Calc[i].atom->number,atoms[FA->num_atm[VC->Calc[i].atom->number]].name,VC->Calc[i].residue->number);
		
		
#if DEBUG_LEVEL > 0
		double com_atm=0.0;
		double Ewall_atm=0.0;
#endif
		
		double radA  = (double)VC->Calc[i].atom->radius;
		double radoA = radA + Rw;
		
		double SAS = 4.0*PI*radoA*radoA;
   		double surfA = SAS;
		
		if(FA->useacs && atoms[atomzero].acs < 0.0){
			// accessible contact surface with solvent/atom		
			// ACS = Total surface area - surface areas of bonded contacts (atoms with a bond/angle between them)
			atoms[atomzero].acs = surfA;
		
			int currindex = VC->ca_index[i];
			
			while(currindex != -1) {
				int atomcont = FA->num_atm[VC->Calc[VC->ca_rec[currindex].atom].atom->number];
				int intramolecular = atoms[atomcont].ofres == atoms[atomzero].ofres;
				
				// get first atom of residue
				int fatm = residue[rnum].fatm[0];
				
				if(intramolecular && residue[rnum].bonded[atomcont-fatm][atomzero-fatm] >= 0)
				{
					atoms[atomzero].acs -= VC->ca_rec[currindex].area;
				}
				
				currindex = VC->ca_rec[currindex].prev;
			}
			
			if(atoms[atomzero].acs < 0.0){ atoms[atomzero].acs = 0.0f; }
			//printf("after ACS=%.3f\n", ACS);
		}
		
#if DEBUG_LEVEL > 0
		cfs_atom.sas = 0.0;
#endif

		if(atoms[atomzero].ncons > 0){

			for(int j=0;j<atoms[atomzero].ncons;j++){
	
				/*
				  double radC = atoms[atomzero].number==atoms[atomzero].cons[j]->anum1 ? 
				  (double)atoms[FA->num_atm[atoms[atomzero].cons[j]->anum2]].radius:
				  (double)atoms[FA->num_atm[atoms[atomzero].cons[j]->anum1]].radius;
				*/
				
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
		printf("ATOM :: RES C RNUM  ANUM  T  RAD ::   COMPL  (W)  DIST   AREA ::     CF.COM     CF.WAL\n");
		printf("----------------------------------------------------------------------------------\n");
		printf("ATOM :: %3s %c %4d %5d %2d %4.2f\n",
		       residue[atoms[atomzero].ofres].name,
		       residue[atoms[atomzero].ofres].chn,
		       residue[atoms[atomzero].ofres].number,
		       atoms[atomzero].number,
		       atoms[atomzero].type,
		       atoms[atomzero].radius);
#endif

		int contnum = 0;  // number of contacts (excluding bloops away atoms)
		int currindex = VC->ca_index[i];
		
		while(currindex != -1) {
			
			double radB  = (double)VC->Calc[VC->ca_rec[currindex].atom].atom->radius;
			// double radoB = radB + Rw;
			// double surfB = 4.0*PI*radoB*radoB;
			
			double rAB   = radA+radB;
			
			//double complementarity = 0.0;
			double area = VC->ca_rec[currindex].area;
			
			struct energy_matrix* energy_matrix = &FA->energy_matrix[(VC->Calc[i].atom->type-1)*FA->ntypes +
										 (VC->Calc[VC->ca_rec[currindex].atom].atom->type-1)];

			//double yval = get_yval(energy_matrix,area/((surfA+surfB)/2.0));
			// always use normalized areas in density functions
			double yval = get_yval(energy_matrix,area/surfA);
			
			// initializes lig and target sol_val to 0.0
            // these values will be modified if the contact
            // implies an atom of the ligand and target
            double lig_sol_val = 0.0, target_sol_val = 0.0;

			
			// the following if/else block is used to get proper ligand's contact with water in lig_sol_val
			// and the proper target's contact with water in tar_sol_val
			// Ideally, there would be only one of the two atoms in contact that is from the ligand

            if(!strcmp(VC->Calc[i].residue->name, FA->resligand->name) 	&&
               VC->Calc[i].residue->chn == FA->resligand->chn 			&&
               VC->Calc[i].residue->number == FA->resligand->number  	)
			{
				energy_matrix = &FA->energy_matrix[(VC->Calc[i].atom->type-1)*FA->ntypes + (FA->ntypes-1)];
				lig_sol_val = get_yval(energy_matrix,area/surfA);

				energy_matrix = &FA->energy_matrix[(VC->Calc[VC->ca_rec[currindex].atom].atom->type-1)*FA->ntypes + (FA->ntypes-1)];
				target_sol_val = get_yval(energy_matrix,area/surfA);
			}
            else if(!strcmp(VC->Calc[VC->ca_rec[currindex].atom].residue->name, FA->resligand->name)	&&
                    VC->Calc[VC->ca_rec[currindex].atom].residue->chn == FA->resligand->chn 			&&
                    VC->Calc[VC->ca_rec[currindex].atom].residue->number == FA->resligand->number   	)
			{
				energy_matrix = &FA->energy_matrix[(VC->Calc[VC->ca_rec[currindex].atom].atom->type-1)*FA->ntypes + (FA->ntypes-1)];
				lig_sol_val = get_yval(energy_matrix,area/surfA);

				energy_matrix = &FA->energy_matrix[(VC->Calc[i].atom->type-1)*FA->ntypes + (FA->ntypes-1)];
				target_sol_val = get_yval(energy_matrix,area/surfA);
			}


			SAS -= area;
			
			// number of contacts counter
			contnum++;
			

#if DEBUG_LEVEL > 0
			cfs_atom.com  =  0.0;
			cfs_atom.wal  =  0.0;
#endif
			// atom in contact with atom zero
			int atomcont = FA->num_atm[VC->Calc[VC->ca_rec[currindex].atom].atom->number];
			
			int intramolecular = 0;
			int intraresidue = 0;
			if(atoms[atomcont].ofres == atoms[atomzero].ofres){
				intraresidue = 1;
				intramolecular = 1;
				
			}else if(residue[atoms[atomcont].ofres].type == 0 &&
				 residue[atoms[atomzero].ofres].type == 0){

				intramolecular = 1;
			}
			
			// get first atom of residue
			int fatm = residue[rnum].fatm[0];
		
			// is contact atom bonded to atom zero
			// if YES, skip contact atom
			if(intraresidue)
			{
				// always skip atoms forming a bond or angle with each other
				if(residue[rnum].bonded[atomcont-fatm][atomzero-fatm] >= 0)
				{

#if DEBUG_LEVEL > 2
					printf("    (B) %3s %c %4d %5d %2d %4.2f :: %7.4f (%s) %6.2f %6.2f :: %10.3f %10.3f\n",
					       VC->Calc[VC->ca_rec[currindex].atom].residue->name,
					       VC->Calc[VC->ca_rec[currindex].atom].residue->chn,
					       VC->Calc[VC->ca_rec[currindex].atom].residue->number,
					       VC->Calc[VC->ca_rec[currindex].atom].atom->number,
					       VC->Calc[VC->ca_rec[currindex].atom].atom->type,
					       VC->Calc[VC->ca_rec[currindex].atom].atom->radius,
					       
					       yval, energy_matrix->weight ? "Y": "N",
					       VC->ca_rec[currindex].dist,
					       VC->ca_rec[currindex].area,
					       
					       cfs_atom.com, cfs_atom.wal);

#endif
					currindex = VC->ca_rec[currindex].prev;
					continue;	  
				}
			}
			
			if(FA->contacts[VC->Calc[VC->ca_rec[currindex].atom].atom->number]){
				//printf("%d already calculated\n",VC->Calc[VC->ca_rec[currindex].atom].atom->number );
				currindex = VC->ca_rec[currindex].prev;
				continue;
			}
			
			// covalently bonded flag
			bool covalent = false;
			
			// number of constraints for contact atom
			int nconscont = atoms[atomcont].ncons;
      			double dist_opt = 0.0;
			
			// do contacting atoms have the same constraint
			constraint* cons = NULL;
			if(nconszero > 0 && nconscont > 0){
				for(int j=0;j<nconszero;j++){
					for(int k=0;k<nconscont;k++){
						if(atoms[atomzero].cons[j]->id == atoms[atomcont].cons[k]->id){
							cons = &FA->constraints[atoms[atomzero].cons[j]->id];
							break;
						}
					}
					if(cons != NULL){break;}
				}
				
				if(cons != NULL){
					
					if(cons->type == 1){
						
						covalent = true;
						dist_opt = cons->bond_len;
						
						cfs->con -= KDIST * GetValueFromGaussian(VC->ca_rec[currindex].dist,dist_opt,cons->max_dist);						
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
				
				double Ewall = KWALL*(pow(VC->ca_rec[currindex].dist,-12.0)-pow(permea*rAB,-12.0));
				
				cfs->wal += Ewall;

#if DEBUG_LEVEL > 0
				Ewall_atm += Ewall_atm;
				cfs_atom.wal += Ewall;
#endif
				// ligand intramolecular clash exceeds threshold
				// add an entry in the dee elimination
				if(intramolecular && type == 1 && Ewall > DEE_WALL_THRESHOLD){
					intraclashes.push_back(std::pair<int,int>(atomzero-fatm,atomcont-fatm));
				}
				
				// Treat everything as rigid
				if ( VC->ca_rec[currindex].dist <= dee_clash*rAB ) { cfs->rclash=1; }

			}
			
			
			if( !covalent ){
				
				if(FA->intramolecular || !intramolecular) {
					
					double contribution = 0.0;
					double ligsolcontribution = 0.0, tarsolcontribution = 0.0;
					if(energy_matrix->weight){
						if(FA->normalize_area){
							contribution = yval*area/surfA;
							ligsolcontribution = lig_sol_val*area/surfA;
							tarsolcontribution = target_sol_val*area/surfA;
						}else{
							contribution = yval*area;
							ligsolcontribution = lig_sol_val*area;
							tarsolcontribution = target_sol_val*area;
						}
					}else{
						contribution = yval;
						ligsolcontribution = lig_sol_val;
						tarsolcontribution = target_sol_val;
					}
					
					if(FA->useacs){
						//printf("USE ACS\n");
						//printf("default contribution=%.3f\n", contribution);
						contribution = contribution * atoms[atomzero].acs/surfA * FA->acsweight;
						ligsolcontribution = ligsolcontribution * atoms[atomzero].acs/surfA * FA->acsweight;
						tarsolcontribution = tarsolcontribution * atoms[atomzero].acs/surfA * FA->acsweight;
						//printf("after contribution=%.3f\n", contribution);
					}
					
					cfs->com += contribution;
					cfs->ligsolv += ligsolcontribution;
					cfs->tarsolv += tarsolcontribution;
					
#if DEBUG_LEVEL > 0
					cfs_atom.com += contribution;
#endif
					
					FA->contributions[(VC->Calc[i].atom->type-1)*FA->ntypes+(VC->Calc[VC->ca_rec[currindex].atom].atom->type-1)] += contribution;
					if((VC->Calc[i].atom->type-1) != (VC->Calc[VC->ca_rec[currindex].atom].atom->type-1))
						FA->contributions[(VC->Calc[VC->ca_rec[currindex].atom].atom->type-1)*FA->ntypes+(VC->Calc[i].atom->type-1)] += contribution;
					
				}
				/*
				  else{
				  printf("skipped intramolecular contact: %d %d\n",
				  atoms[atomzero].number, atoms[atomcont].number); //VC->Calc[VC->ca_rec[currindex].atom].atom->number);
				  } 
				*/
			}
			
			/*
			// generate surface area matrix
			matrix[VC->Calc[i].atom->type][VC->Calc[VC->ca_rec[currindex].atom].atom->type] += area;
			
			if(VC->Calc[i].atom->type != VC->Calc[VC->ca_rec[currindex].atom].atom->type){
			matrix[VC->Calc[VC->ca_rec[currindex].atom].atom->type][VC->Calc[i].atom->type] += area;
			}
			*/
					
#if DEBUG_LEVEL > 0
			printf("        %3s %c %4d %5d %2d %4.2f :: %7.4f (%s) %6.2f %6.2f :: %10.3f %10.3f\n",
			       VC->Calc[VC->ca_rec[currindex].atom].residue->name,
			       VC->Calc[VC->ca_rec[currindex].atom].residue->chn,
			       VC->Calc[VC->ca_rec[currindex].atom].residue->number,
			       VC->Calc[VC->ca_rec[currindex].atom].atom->number,
			       VC->Calc[VC->ca_rec[currindex].atom].atom->type,
			       VC->Calc[VC->ca_rec[currindex].atom].atom->radius,
						
			       yval, energy_matrix->weight ? "Y": "N",
			       VC->ca_rec[currindex].dist,
			       VC->ca_rec[currindex].area,
			       
			       cfs_atom.com, cfs_atom.wal);
			
#endif

			// skip to next contact
			currindex = VC->ca_rec[currindex].prev;
		}
		
		//    printf("Atom[%d]=%d has %d contacts\n",VC->Calc[i].,VC->Calc[i].atom->number,contnum);
		//    printf("Atom[%d] COM=[%8.2f]\tWAL=[%8.2f]\n",VC->Calc[i].atomnum,com_atm,Ewall_atm);
					
		if(SAS < 0.0){ SAS = 0.0; }
		cfs->totsas += SAS;
		
		double contribution = 0.0;
		if(FA->solventterm){
			contribution = (double)FA->solventterm * SAS;
			//printf("SP: multiply ST=%.3f with SAS.area=%.3f\n", (double)FA->solventterm, SAS);
		} else {
			struct energy_matrix* energy_matrix = &FA->energy_matrix[(VC->Calc[i].atom->type-1)*FA->ntypes +
										 (FA->ntypes-1)];
			//printf("type1: %d\ttype2: %d\n", energy_matrix->type1, energy_matrix->type2);
			
			double yval = get_yval(energy_matrix,SAS/surfA);
			
			if(energy_matrix->weight){
				if(FA->normalize_area){
					contribution = yval * SAS / surfA;
				}else{
					contribution = yval * SAS;
				}
				//printf("Weight: multiply yval=%.3f by SAS.area=%.3f\n", yval, SAS);
			}
			else {
				contribution = yval;
				/*
				  if(VC->Calc[i].type == 3){
				  printf("Density: add yval=%.3f for norm.SAS=%.3f for atom %d\n",
				  yval, SAS/surfA, VC->Calc[i].atomnum);
				  }
				*/
			}
		}
		
		if(FA->useacs){
			contribution = contribution * atoms[atomzero].acs/surfA * FA->acsweight;
		}
		
		cfs->sas += contribution;
		
		FA->contributions[(VC->Calc[i].atom->type-1)*FA->ntypes + (FA->ntypes-1)] += contribution;
		FA->contributions[(FA->ntypes-1)*FA->ntypes + (VC->Calc[i].atom->type-1)] += contribution;
		
		FA->contacts[VC->Calc[i].atom->number] = 1;
		
#if DEBUG_LEVEL > 1
		printf("CF.SAS is %.3f for %d contacts with contribution %.3f\n", 
		       SAS, contnum, contribution); 
#endif
		
	}
	

	// Penalize Freesurf.
  
	//printf("FREESURF(SAS)=[%8.2f]\n",SAStot);
	//printf("FINAL SUM COM=[%8.2f]\tWAL=[%8.2f]\n",com,Ewall);
  
	//print_surfmat(matrix,"surf.mat");


/*
  #if DEBUG_LEVEL > 0
  printf("\n");
  printf("CF.sum = %.3f\n", cfs->com + cfs->sas + cfs->wal);
  printf("CF.com = %.3f\n", cfs->com);
  printf("CF.sas = %.3f\n", cfs->sas);
  printf("CF.wal = %.3f\n", cfs->wal);
  getchar(); 
  #endif
*/
	//getchar();
	
	for(int i=0;i<FA->atm_cnt_real;i++){
		if(VC->Calc[i].score){
			VC->Calc[i].atom = NULL;
		}
	}

	if(!FA->vindex){ free(VC->box); }
	
	return(0.0);
  
}

double get_yval(struct energy_matrix* energy_matrix, double relative_area)
{
	double yval = 0.0;
	
	// a single value in matrix (weighted by area)
	if(energy_matrix->weight)
		yval = energy_matrix->energy_values->y;
	else { // density function
		struct energy_values* xyval = energy_matrix->energy_values;

		while(xyval->next_value != NULL && relative_area > xyval->next_value->x){
			/*
			  printf("x=%.3f next_value.x=%.3f next_value.y=%.3f\n",
			  xyval->x, xyval->next_value->x, xyval->next_value->y);
			*/
			xyval = xyval->next_value;
		}
		
		if(xyval->x > relative_area){
			// no left bound data
			yval = 0.0;
		}else if(xyval->next_value == NULL){
			// no right bound data
			yval = xyval->y;
		}else{
			yval = xyval->y + 
				( relative_area - xyval->x ) / (xyval->next_value->x - xyval->x ) *
				( xyval->next_value->y - xyval->y );
		}
		
		/*
		  if(energy_matrix->type2 == 40){
		  printf("stopped at x=%.3f with y=%.3f\n", xyval->x, xyval->y);
		  if(xyval->next_value != NULL){
		  printf("next is x=%.3f with y=%.3f\n", xyval->next_value->x, xyval->next_value->y);
		  }
		  printf("prob func. yval=%.3f for relative_area %.3f for [%d][%d]\n", yval, relative_area,
		  energy_matrix->type1, energy_matrix->type2);
		  printf("calculated y=%.3f\n", yval);
		  getchar();
		  }
		*/
	}

	return yval;
}
