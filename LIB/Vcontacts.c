#include "Vcontacts.h"
#include "boinc.h"

// Vcontacts calculates the SAS only for the residue sent in argument
int Vcontacts(FA_Global* FA,atom* atoms,resid* residue,VC_Global* VC)
{  
    
	//VC->planedef = 'X';  // extended radical plane (default)
	VC->planedef = 'R';  // radical plane
	//VC->planedef = 'B';  // bisection
    
	// initialize contact atom index
	VC->ca_recsize = 5*FA->atm_cnt_real;
	VC->ca_rec = (ca_struct*)malloc(VC->ca_recsize*sizeof(ca_struct));
	//printf("allocating %p\n",VC->ca_rec);
    
	if(!VC->ca_rec) {
		fprintf(stderr,"ERROR: memory allocation error for ca_rec\n"); 
		Terminate(2);
	}
    
	for(int i=0; i<FA->atm_cnt_real; ++i) {
		VC->Calc[i].vol = 0.0;
		VC->ca_index[i] = -1;   //initialize pointer array
		VC->seed[i*3] = -1;
	}
    
	// protein atoms to boxes in cubic grid
	index_protein(FA,atoms,residue,VC->Calc,&VC->box,VC->Calclist,&VC->dim,FA->atm_cnt_real);
    
	VC->numcarec=0;
	// calc volumes for all protein atoms
    
	return calc_region(FA,VC,atoms,FA->atm_cnt_real);
    
}

/***************************
 * subroutine calc_region
 ***************************/

// this subroutine calculates the contact areas (SAS) for a given set of atoms.
// Here, the set of atoms is the entire protein.
// needs global variable 'dim'.

int calc_region(FA_Global* FA,VC_Global* VC,atom* atoms,int atmcnt)
{
	int    i;        // atom counter
	int    atomzero; // current center atom
	int    boxi;
	int    NC;       // number of contacts around atomzero
	int    NV;       // number of vertices in polyhedron around atomzero
	float  rado;     // radius of atomzero PLUS radius of water
	char   surfatom; // atom type, 'I' internal, 'S' surface
    
	/*
     printf("================================\n");
     printf("Calculating SAS for residue [%d]\n", resnum);
     */
    
	for(i=0;i<atmcnt;++i) {
		// ============= atom contact calculations =============
		atomzero = VC->Calclist[i];
		boxi = VC->Calc[atomzero].boxnum;
        
		if(!VC->Calc[atomzero].score){continue;}
        
		//printf("Get_contacts for %d\n",VC->Calc[atomzero].atomnum);   
		rado = VC->Calc[atomzero].radius + Rw;
        
		NC = get_contlist4(atoms,atomzero, VC->contlist, atmcnt, rado, VC->dim, VC->Calc, VC->Calclist, VC->box,VC->ca_rec, VC->ca_index);
        
		// invalid write/read when NC = 0
		// because planeA is negative subscript
		if(NC == 0){
			VC->Calc[atomzero].SAS = -1.0;
			continue;
		}
        
		NV = voronoi_poly2(VC,atomzero, VC->cont, rado, NC, VC->contlist); 
        
		// could not generate polyhedron
		// because of clashing atom
		if(NV == -1){
			return(-1);
		}
        
		surfatom = order_faces(atomzero, VC->poly, VC->centerpt, rado, NC, NV, VC->cont, VC->ptorder);    
        
		calc_areas(VC->poly, VC->centerpt, rado, NC, NV, VC->cont, VC->ptorder, &VC->Calc[atomzero]);                
        
		save_areas(VC->cont, VC->contlist, NC, atomzero, VC->Calc,&VC->ca_recsize, &VC->numcarec, &VC->ca_rec, VC->ca_index);
        
		min_areas(VC->ca_rec, VC->Calc, &VC->Calc[atomzero], FA->vcontacts_self_consistency);
        
	}
    
    
	return(0);
}


/******************************
 * subroutine voronoi_poly2
 * created 08/07/2001  BJM
 ******************************/

int voronoi_poly2(VC_Global *VC,int atomzero, plane cont[], float rado, 
                  int NC, const contactlist* contlist)
{
	int    cai;          // contact atom counter
	atomsas *ca_ptr; // pointer to pdb atom
	double atomdist;     // distance to atom
	double planedist;    // distance to plane
	double mindist;      // distance to closest plane
	int    planeA;       // closest plane to origin
	int    planeB;       // second plane, with planeA defines closest edge
	int    planeC;       // new intersection plane for edge (endpt)
	int    oldplaneC;    // old intersection plane for edge (startpt)
	double vt;           // vector parameter 't' for line x=x'+lt, y=y'+mt, z=z'+nt
	double vtmin;        // minimum value for vt 
	double vtdiv;        // check for division by zero in vt calculation 
	double temppt[3];
	vertex *vstp; // pointer to start vertex, coordinates
	double *V;          // pointer to edge vector
	int    startedge;
	int    edgenum;
	int    vn = 0;
	char   edgeflag;
	int    edgei;      // edge counter
	int    vi, vj;     // vertices counters
	double arcpt0[3], arcpt1[3];
	int    testpA, testpB;
	double testvalA, testvalB;
	char   arcflag = 'N';
    
	// failsafe variables:
	char   recalc;       // flag if hull is being recalculated (orig. unbounded)
	float  origcoor[3];  // original pdb coordinates for atom. 
    
	recalc = 'N';
RESTART:
	planeA = -1;
	planeB = -1;
	planeC = -1;
    
	origcoor[0] = 0.0f;
	origcoor[1] = 0.0f;
	origcoor[2] = 0.0f;
    
	/* generate planes of contact with D = planedist */
	mindist = 9.9e+9;
	for(cai=0; cai<NC; ++cai) {
		ca_ptr = &VC->Calc[contlist[cai].index];
		atomdist = contlist[cai].dist;
        
		if(VC->planedef == 'B') {  // bisection - original Voronoi procedure
			planedist = atomdist/2.0;
		} else if(VC->planedef == 'R') { // radical plane (Gellatly and Finney) - default.
			planedist = (atomdist*atomdist + (rado-Rw)*(rado-Rw) - (ca_ptr->radius)*(ca_ptr->radius))/(2*atomdist);
		} else { // extended radical plane (McConkey et al). 
			planedist = (atomdist*atomdist + rado*rado - (Rw + ca_ptr->radius)*(Rw + ca_ptr->radius))/(2*atomdist);
		}
        
		cont[cai].Ai[0] = (ca_ptr->coor[0] - VC->Calc[atomzero].coor[0])/atomdist;
		cont[cai].Ai[1] = (ca_ptr->coor[1] - VC->Calc[atomzero].coor[1])/atomdist;
		cont[cai].Ai[2] = (ca_ptr->coor[2] - VC->Calc[atomzero].coor[2])/atomdist;
		cont[cai].Ai[3] = -planedist;
		cont[cai].dist  = fabs(planedist);
		cont[cai].index = contlist[cai].index;
		cont[cai].flag = 'X'; // initialize contact flags to 'no contact' 
        
		//printf("plane[%d] with dist[%8.3f] for atom %d\n",cai,cont[cai].dist,VC->Calc[cont[cai].index].atomnum);
		// set plane0 as closest plane
		if(cont[cai].dist < mindist) {
			mindist = cont[cai].dist;
			planeA = cai;
			//printf("min plane=%d\tmin dist=%8.3f for atom %d\n",planeA,cont[cai].dist,VC->Calc[cont[cai].index].atomnum);
		}
	}
    
	// add four planes surrounding atom, outer limit for voronoi polyhedron
	cont[NC].Ai[0] = 0.707;
	cont[NC].Ai[1] = 1.0;
	cont[NC].Ai[2] = 0.0;
	cont[NC].Ai[3] = -10.0;
    
	cont[NC+1].Ai[0] = 0.707;
	cont[NC+1].Ai[1] = -1.0;
	cont[NC+1].Ai[2] = 0.0;
	cont[NC+1].Ai[3] = -10.0;
    
	cont[NC+2].Ai[0] = -0.707;
	cont[NC+2].Ai[1] = 0.0;
	cont[NC+2].Ai[2] = 1.0;
	cont[NC+2].Ai[3] = -10.0;
    
	cont[NC+3].Ai[0] = -0.707;
	cont[NC+3].Ai[1] = 0.0;
	cont[NC+3].Ai[2] = -1.0;
	cont[NC+3].Ai[3] = -10.0;
    
	// get starting vertex from seed or calc new vertex
    
	get_firstvert(VC->seed, cont, &planeA, &planeB, &planeC, NC, atomzero);
    
	solve_3x3(cont[planeA].Ai, cont[planeB].Ai, cont[planeC].Ai, temppt);
    
	// add first vertex to vertex list
	add_vertex(VC->poly, 0, temppt, planeA, planeB, planeC);
    
	// flag contacts as present
	cont[planeA].flag = 'Y';
	cont[planeB].flag = 'Y';
	cont[planeC].flag = 'Y';
    
    
	// calculate edge vectors
	add_vedge(VC->vedge, 0, cont, planeA, planeB, planeC, VC->poly, 0);
	add_vedge(VC->vedge, 1, cont, planeB, planeC, planeA, VC->poly, 0);
	add_vedge(VC->vedge, 2, cont, planeC, planeA, planeB, VC->poly, 0);
    
	startedge = 0;
	edgenum = 3;
	vn = 1;
    
	/* --------------------------------------------------- */
	/* Generate new polyhedron points from edge vectors    */
	/* --------------------------------------------------- */
    
	while(1) {
		// get next unfinished vector = startedge
		while((VC->vedge[startedge].endpt >= 0) && ((edgenum-startedge) > 0)) {
			++startedge;
		}
		if((edgenum-startedge) <= 0) {
			// all edges are done, polyhedron complete.
			break;
		}
        
		vtmin = 9.9e+9; // dummy value
		vstp = &VC->poly[VC->vedge[startedge].startpt];
		V = VC->vedge[startedge].V;
		planeA = VC->vedge[startedge].plane[0];
		planeB = VC->vedge[startedge].plane[1];
		oldplaneC = VC->vedge[startedge].startplane;
		planeC = -1;
        
		// get closest positive intersection point
		for(cai=0; cai<NC; ++cai) {
			// check if contact is to be done - for now, do all.
			if((cai != planeA) && (cai != planeB) && (cai != oldplaneC)) {
				vtdiv = (cont[cai].Ai[0]*V[0] +cont[cai].Ai[1]*V[1] +cont[cai].Ai[2]*V[2]);
				if(vtdiv != 0.0) {
					vt = -(cont[cai].Ai[0]*vstp->xi[0] +cont[cai].Ai[1]*vstp->xi[1] 
					       +cont[cai].Ai[2]*vstp->xi[2] +cont[cai].Ai[3])/vtdiv;
					if((vt < vtmin) && (vt > 0)) {
						vtmin = vt;
						planeC = cai;
					}
				}
			}
		}
		VC->poly[vn].xi[0] = vstp->xi[0] + vtmin*V[0];
		VC->poly[vn].xi[1] = vstp->xi[1] + vtmin*V[1];
		VC->poly[vn].xi[2] = vstp->xi[2] + vtmin*V[2];
        
		// if point is outside sphere, check vs. external planes
		if((VC->poly[vn].xi[0]*VC->poly[vn].xi[0] + VC->poly[vn].xi[1]*VC->poly[vn].xi[1] 
		    + VC->poly[vn].xi[2]*VC->poly[vn].xi[2]) > rado*rado) {
			for(cai=NC; cai<NC+4; ++cai) {
				// check if contact is to be done - for now, do all.
				if((cai != planeA) && (cai != planeB) && (cai != oldplaneC)) {
					vtdiv = (cont[cai].Ai[0]*V[0] +cont[cai].Ai[1]*V[1] +cont[cai].Ai[2]*V[2]);
					if(vtdiv != 0.0) {
						vt = -(cont[cai].Ai[0]*vstp->xi[0] +cont[cai].Ai[1]*vstp->xi[1] 
						       +cont[cai].Ai[2]*vstp->xi[2] +cont[cai].Ai[3])/vtdiv;
						if((vt < vtmin) && (vt > 0)) {
							vtmin = vt;
							planeC = cai;
						}
					}
				}
			}
			VC->poly[vn].xi[0] = vstp->xi[0] + vtmin*V[0];
			VC->poly[vn].xi[1] = vstp->xi[1] + vtmin*V[1];
			VC->poly[vn].xi[2] = vstp->xi[2] + vtmin*V[2];
		}
        
		add_vertex(VC->poly, vn, VC->poly[vn].xi, planeA, planeB, planeC);
		VC->vedge[startedge].endpt = vn;
		VC->vedge[startedge].endplane = planeC;
        
		//flag contact as present
		cont[planeC].flag = 'Y';
        
		// ========  ADD EDGES  ========
        
		// check edge (planeA, planeC)
		edgeflag = 'Y';
		edgei = startedge+1;
		while(edgei < edgenum) {
			if(((VC->vedge[edgei].plane[0] == planeA)&&(VC->vedge[edgei].plane[1] == planeC)) ||
			   ((VC->vedge[edgei].plane[0] == planeC)&&(VC->vedge[edgei].plane[1] == planeA))) {
				// already on list, add current vertex as endpt
				VC->vedge[edgei].endpt = vn;
				VC->vedge[edgei].endplane = planeB;
				edgeflag = 'N';
				break;
			}
			++edgei;
		}
		if(edgeflag == 'Y') { // add edge
			add_vedge(VC->vedge, edgenum, cont, planeA, planeC, planeB, VC->poly, vn);
			++edgenum;
		}
        
		// check edge (planeB, planeC)
		edgeflag = 'Y';
		edgei = startedge+1;
		while(edgei < edgenum) {
			if(((VC->vedge[edgei].plane[0] == planeB)&&(VC->vedge[edgei].plane[1] == planeC)) ||
			   ((VC->vedge[edgei].plane[0] == planeC)&&(VC->vedge[edgei].plane[1] == planeB))) {
				// already on list, add current vertex as endpt
				VC->vedge[edgei].endpt = vn;
				VC->vedge[edgei].endplane = planeA;
				edgeflag = 'N';
				break;
			}
			++edgei;
		}
		if(edgeflag == 'Y') { // add edge
			add_vedge(VC->vedge, edgenum, cont, planeB, planeC, planeA, VC->poly, vn);
			++edgenum;
            
			// ===== failsafe - if solution is not converging, perturb atom  =====
			// ===== coordinates and recalculate.                            =====
			if(edgenum >= 200) {
				//printf("********* invalid solution for hull, recalculating *********\n");
				VC->seed[atomzero*3] = -1;  // reset to no seed vertex
				origcoor[0] = VC->Calc[atomzero].coor[0];
				origcoor[1] = VC->Calc[atomzero].coor[1];
				origcoor[2] = VC->Calc[atomzero].coor[2];
                
				// perturb atom coordinates
				VC->Calc[atomzero].coor[0] += 0.005f*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
				VC->Calc[atomzero].coor[1] += 0.005f*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
				VC->Calc[atomzero].coor[2] += 0.005f*(float)(2*rand()-RAND_MAX)/(float)RAND_MAX;
                
				// *** NEW ***
                
				// Do not recalc
				if (VC->first) {
					recalc = 'Y';
                    
					// EXCEPT REFERENCE SOLUTION (FIRST CALL TO VCT)
					// Never recalculate because solution that do not converge are clashing solutions
					// Those individuals would not survive in evolution
					goto RESTART;
				}
                
				// Abort immediately Scoring
                
				return -1;
			}
		}
		++vn;
	}
    
	/*--------------------------------------------------*/
	/*  now have voronoi polyhedron around given atom.  */
	/*  remove vertices outside of sphere, and          */
	/*  calculate intersection points with sphere.      */ 
	/*--------------------------------------------------*/
    
	// flag edges that may cross sphere boundary
	for(edgei=0; edgei<edgenum; ++edgei) {
		if((rado < VC->poly[VC->vedge[edgei].startpt].dist) || (rado < VC->poly[VC->vedge[edgei].endpt].dist)) {
			// one or both vertices fall outside of sphere
			arcflag = 'Y';
			VC->vedge[edgei].arc = '?';
		} else {
			VC->vedge[edgei].arc = 'X';
		}
	}
    
	// calculate new arc points
	for(edgei=0; edgei<edgenum; ++edgei) {
		if(VC->vedge[edgei].arc != 'X') {
			if(solve_2xS(&cont[VC->vedge[edgei].plane[0]], &cont[VC->vedge[edgei].plane[1]], rado, arcpt0, arcpt1) == -1) {
				VC->vedge[edgei].arc = 'X';  // mark edge as no associated arc point
				continue;
			}
			// test new arc points vs. adjacent planes, add if ok.
			testpA = VC->vedge[edgei].startplane;
			testpB = VC->vedge[edgei].endplane;
			testvalA = cont[testpA].Ai[0]*arcpt0[0] + cont[testpA].Ai[1]*arcpt0[1]
            + cont[testpA].Ai[2]*arcpt0[2] + cont[testpA].Ai[3];
			testvalB = cont[testpB].Ai[0]*arcpt0[0] + cont[testpB].Ai[1]*arcpt0[1] 
            + cont[testpB].Ai[2]*arcpt0[2] + cont[testpB].Ai[3];
			if((testvalA < 0.0) && (testvalB < 0.0)) { // point is good
				add_vertex(VC->poly, vn, arcpt0, VC->vedge[edgei].plane[0], VC->vedge[edgei].plane[1], -1);
				VC->poly[vn].dist = rado;
				++vn;
			}
			testvalA = cont[testpA].Ai[0]*arcpt1[0] + cont[testpA].Ai[1]*arcpt1[1] 
            + cont[testpA].Ai[2]*arcpt1[2] + cont[testpA].Ai[3];
			testvalB = cont[testpB].Ai[0]*arcpt1[0] + cont[testpB].Ai[1]*arcpt1[1] 
            + cont[testpB].Ai[2]*arcpt1[2] + cont[testpB].Ai[3];
			if((testvalA < 0.0) && (testvalB < 0.0)) { // point is good
				add_vertex(VC->poly, vn, arcpt1, VC->vedge[edgei].plane[0], VC->vedge[edgei].plane[1], -1);
				VC->poly[vn].dist = rado;
				++vn;
			}
		}
	}
    
	// reduce poly vertex list
	vj=0;
	for(vi=0; vi<vn; ++vi) {
		if(VC->poly[vi].dist <= rado) {
			VC->poly[vj] = VC->poly[vi];
			++vj;
		}
	}
	vn = vj;
    
	// ----- calculate center points and mark engulfing contacts -----
	for(cai=0; cai<NC; ++cai) {
		VC->centerpt[cai].xi[0] = -cont[cai].Ai[0]*cont[cai].Ai[3];
		VC->centerpt[cai].xi[1] = -cont[cai].Ai[1]*cont[cai].Ai[3];
		VC->centerpt[cai].xi[2] = -cont[cai].Ai[2]*cont[cai].Ai[3];
		VC->centerpt[cai].dist = fabs(cont[cai].Ai[3]);
		if(cont[cai].Ai[3] > 0.0) {
			cont[cai].flag = 'E';
		}
	}
    
	if(recalc == 'N') {
		save_seeds(VC->seed,cont, VC->poly, vn, atomzero);
	} else {
		// reset atom coordinates to original values
		VC->Calc[atomzero].coor[0] = origcoor[0];
		VC->Calc[atomzero].coor[1] = origcoor[1];
		VC->Calc[atomzero].coor[2] = origcoor[2];
	}
    
	return(vn);
}



/****************************
 * subroutine get_firstvert
 ****************************/

void get_firstvert(const int* seed,const plane* cont, int *planeA, int *planeB, int *planeC, int NC, int atomzero)
{
	int seedi;
	int cai;
	double mindist;
	double ptA[3];
	double vectA[4];     // 4 values so it can be used as a plane equation as well
	double vt;           // vector parameter 't' for line x=x'+lt, y=y'+mt, z=z'+nt
	double vtmin;        // minimum value for vt 
	double vtdiv;        // check for division by zero in vt calculation 
	double temppt[3];
	double ptdist;
    
	// get previous seed vertex if present
	seedi = atomzero*3;
	*planeA = -1;
	*planeB = -1;
	*planeC = -1;
	if(seed[seedi] != -1) {
		for(cai=0; cai<NC; ++cai) {
			if(cont[cai].index == seed[seedi]) {
				*planeA = cai;
			} else if(cont[cai].index == seed[seedi+1]) {
				*planeB = cai;
			} else if(cont[cai].index == seed[seedi+2]) {
				*planeC = cai;
			}
		}
	}
    
	if((*planeA != -1)&&(*planeB != -1)&&(*planeC != -1)) {
		return;
	} else {
        
		// -------------  find initial edge, on plane closest to origin  -------------
		mindist = 9.9e+9;  // dummy value
		for(cai=0; cai<NC; ++cai) {
			if(cont[cai].dist < mindist) {
				mindist = cont[cai].dist;
				*planeA = cai;
			}
		}
        
		ptA[0] = 0.0;
		ptA[1] = 0.0;
		ptA[2] = 0.0;
		mindist = 9.9e+9;  // dummy value
        
		for(cai=0; cai<NC+4; ++cai) {
			if(cai != *planeA) {
				/*
                 printf("*planeA: %d\tcai: %d\tNC: %d\n", *planeA, cai, NC);
                 if(*planeA == -1) { getchar(); }
                 */
				vectA[0] = cont[*planeA].Ai[1]*cont[cai].Ai[2] - cont[*planeA].Ai[2]*cont[cai].Ai[1];
				vectA[1] = cont[*planeA].Ai[2]*cont[cai].Ai[0] - cont[*planeA].Ai[0]*cont[cai].Ai[2];
				vectA[2] = cont[*planeA].Ai[0]*cont[cai].Ai[1] - cont[*planeA].Ai[1]*cont[cai].Ai[0];
				vectA[3] = 0.0;
				if(solve_3x3(cont[*planeA].Ai, cont[cai].Ai, vectA, temppt) != -1) {
					ptdist = sqrt(temppt[0]*temppt[0] + temppt[1]*temppt[1] + temppt[2]*temppt[2]);
					if(ptdist < mindist) {
						*planeB = cai;
						mindist = ptdist;
						ptA[0] = temppt[0];
						ptA[1] = temppt[1];
						ptA[2] = temppt[2];
					}
				}
			}
		}
        
		// recalc vector normal to planes A and B
		vectA[0] = cont[*planeA].Ai[1]*cont[*planeB].Ai[2] - cont[*planeA].Ai[2]*cont[*planeB].Ai[1];
		vectA[1] = cont[*planeA].Ai[2]*cont[*planeB].Ai[0] - cont[*planeA].Ai[0]*cont[*planeB].Ai[2];
		vectA[2] = cont[*planeA].Ai[0]*cont[*planeB].Ai[1] - cont[*planeA].Ai[1]*cont[*planeB].Ai[0];
		vectA[3] = 0.0;
        
		// get starting vertex on polyhedron
		vtmin = 9.9e+9;   // dummy value
		for(cai=0; cai<NC+4; ++cai) {
			if((cai != *planeA) && (cai != *planeB)) {
				vtdiv = (cont[cai].Ai[0]*vectA[0] +cont[cai].Ai[1]*vectA[1] +cont[cai].Ai[2]*vectA[2]);
				if(vtdiv != 0.0) {
					vt = -(cont[cai].Ai[0]*ptA[0] +cont[cai].Ai[1]*ptA[1] +cont[cai].Ai[2]*ptA[2] +cont[cai].Ai[3])/vtdiv;
					if(fabs(vt) < vtmin) {
						vtmin = fabs(vt);
						*planeC = cai;
					}
				}
			}
		}
	}
    
	return;
}


/*************************
 * subroutine save_seeds
 *************************/

// saves starting vertices for atoms to be done

void save_seeds(int* seed,const plane* cont, const vertex* poly, int NV, int atomzero)
{
	int vi;
	int seedi;
    
	for(vi=0; vi<NV; ++vi) {
		if(poly[vi].plane[2] != -1) {
			seedi = 3*cont[poly[vi].plane[0]].index;
			if(seed[seedi] == -1) {
				seed[seedi] = atomzero;
				seed[seedi+1] = cont[poly[vi].plane[1]].index;
				seed[seedi+2] = cont[poly[vi].plane[2]].index;
			}
			seedi = 3*cont[poly[vi].plane[1]].index;
			if(seed[seedi] == -1) {
				seed[seedi] = atomzero;
				seed[seedi+1] = cont[poly[vi].plane[0]].index;
				seed[seedi+2] = cont[poly[vi].plane[2]].index;
			}
			seedi = 3*cont[poly[vi].plane[2]].index;
			if(seed[seedi] == -1) {
				seed[seedi] = atomzero;
				seed[seedi+1] = cont[poly[vi].plane[0]].index;
				seed[seedi+2] = cont[poly[vi].plane[1]].index;
			}
		}
	}  
	return;
}



/*******************************
 * subroutine add_vedge
 *******************************/

// adds a new edge to the edge list. Direction of vector is tested so that
// it points away from the startpt, along the body of the polyhedron.

// stores information in variable vedge[edgenum].

void add_vedge(edgevector vedge[], int edgenum, const plane* cont, int plane0, int plane1, 
               int testplane, const vertex* poly, int startpt)
{
	double testpt[3];
	double testval;
    
	vedge[edgenum].V[0] = cont[plane0].Ai[1]*cont[plane1].Ai[2] - cont[plane0].Ai[2]*cont[plane1].Ai[1];
	vedge[edgenum].V[1] = cont[plane0].Ai[2]*cont[plane1].Ai[0] - cont[plane0].Ai[0]*cont[plane1].Ai[2];
	vedge[edgenum].V[2] = cont[plane0].Ai[0]*cont[plane1].Ai[1] - cont[plane0].Ai[1]*cont[plane1].Ai[0];
	vedge[edgenum].startpt = startpt;
	vedge[edgenum].endpt = -1;  // flag, edge not completed.
	vedge[edgenum].plane[0] = plane0;
	vedge[edgenum].plane[1] = plane1;
	vedge[edgenum].startplane = testplane;
	vedge[edgenum].arc = '.'; // dummy value.
	//printf("vedge[edgenum(%d)]\n",edgenum);
    
	// test direction of vector
	testpt[0] = poly[startpt].xi[0] + vedge[edgenum].V[0];
	testpt[1] = poly[startpt].xi[1] + vedge[edgenum].V[1];
	testpt[2] = poly[startpt].xi[2] + vedge[edgenum].V[2];
    
	testval = cont[testplane].Ai[0]*testpt[0] +cont[testplane].Ai[1]*testpt[1] 
    + cont[testplane].Ai[2]*testpt[2] +cont[testplane].Ai[3];
	if(testval > 0.0) { // vector is in wrong direction
		vedge[edgenum].V[0] = -vedge[edgenum].V[0];
		vedge[edgenum].V[1] = -vedge[edgenum].V[1];
		vedge[edgenum].V[2] = -vedge[edgenum].V[2];
	}
	return;
}

/*******************************
 * subroutine add_vertex
 *******************************/

// add polyhedron vertex to local list for current atom contacts

int add_vertex(vertex poly[], int vn, const double* coor, int planeA, int planeB, int planeC)
{ 
	poly[vn].xi[0] = coor[0];
	poly[vn].xi[1] = coor[1];
	poly[vn].xi[2] = coor[2];
	poly[vn].plane[0] = planeA;
	poly[vn].plane[1] = planeB;
	poly[vn].plane[2] = planeC;
	poly[vn].dist = sqrt(coor[0]*coor[0] +coor[1]*coor[1] +coor[2]*coor[2]);
	return(0);
}

/*********************************
 * function order_faces          
 *********************************/

// output is array ptorder, the order of points around each face.
// return values are 'I' for internal atom, 'S' for surface atom.

char order_faces(int atomzero,vertex poly[], const vertex* centerpt, float rado,
                 int NC, int NV, const plane* cont, ptindex ptorder[])
{
	int planeX;              // current plane, 1 to NC
	int planeY;              // second plane, to get adjacent points on polygon
	int surfcount;           // number of points defining a given plane of contact
	int vi, vi2;             // vertices counter
	int tempsi;              // temporary storage for exchanging indices
	double tempcos;          // for exchanging cosines
	double cos10X[MAX_CONT]; // cos of angle between points poly[1],poly[0],and poly[X]
	double temppt[3];        // temp coordinates of new arc point
	char surfatom;           // return value: 'I' internal atom, 'S' surface atom
    
	surfatom = 'I'; // internal atom
	for(vi=0; vi<NV; ++vi) {
		if(poly[vi].plane[2] == -1) {
			surfatom = 'S'; // surface atom 
			break;
		}
	}
    
	//printf("here1\n");
	// for surface calculation only
	// if(surfatom == 'I') return('I');
    
	for(planeX=0; planeX < NC; ++planeX) {
		//printf("planeX=%d\n",planeX);
        
		if(cont[planeX].flag == 'X') { // hidden
			ptorder[planeX].numpts = 0;
			continue;
		}
        
		//printf("here2\n");
		surfcount=0;
		for(vi=0; vi<NV; ++vi) {
			// index all points comprising surface for planeX
			if((poly[vi].plane[0]==planeX) || (poly[vi].plane[1]==planeX) || (poly[vi].plane[2]==planeX)) {
				ptorder[planeX].pt[surfcount] = vi;
				++surfcount;
			}
		}
		//printf("here3\n");
        
		if(surfcount > 2) {
            
			/*-------------------------------------------------------------------*/
			/* three or more surface points, need additional point so no arcs    */
			/* are greater than 180 deg. Make point opposite first pt on an arc. */
			/* (if no points on arc, extra point isn't needed.)                  */
			/*-------------------------------------------------------------------*/
            
			for(vi=0; vi<surfcount; ++vi) {
				if(poly[ptorder[planeX].pt[vi]].plane[2] == -1) { // on arc, calc pt. opposite
                    
					// get coordinates of point
					temppt[0] = 2.0*centerpt[planeX].xi[0] - poly[ptorder[planeX].pt[vi]].xi[0];
					temppt[1] = 2.0*centerpt[planeX].xi[1] - poly[ptorder[planeX].pt[vi]].xi[1];
					temppt[2] = 2.0*centerpt[planeX].xi[2] - poly[ptorder[planeX].pt[vi]].xi[2];
                    
					//keep point if it's valid
					if(test_point(temppt, cont, NC, rado, planeX, -1, -1) == 'Y') {
						poly[NV].xi[0] = temppt[0];
						poly[NV].xi[1] = temppt[1];
						poly[NV].xi[2] = temppt[2];
						poly[NV].plane[0] = planeX;
						poly[NV].plane[1] = -1;
						poly[NV].plane[2] = -1;
						poly[NV].dist = rado;
						ptorder[planeX].pt[surfcount] = NV;
						++surfcount;
						++NV;
					}
					break;
				}
			}
			//printf("here4\n");
            
		}
        
		//printf("here5\n");
        
		ptorder[planeX].numpts = surfcount;
        
		//printf("here6\n");
        
		if(surfcount > 3) {
			// get two points on same line (two common planes).
			// all points already share one plane (planeX), find another.
			if(poly[ptorder[planeX].pt[0]].plane[0] == planeX) {
				planeY = poly[ptorder[planeX].pt[0]].plane[1];
			} else {
				planeY = poly[ptorder[planeX].pt[0]].plane[0];
			}
            
			//find another point on the same line (2 common planes)
			for(vi=1; vi<surfcount; ++vi) {
				if((poly[ptorder[planeX].pt[vi]].plane[0]==planeY) || (poly[ptorder[planeX].pt[vi]].plane[1]==planeY) 
				   || (poly[ptorder[planeX].pt[vi]].plane[2]==planeY)) {
					break;
				}
			}    
            
			//swap index for pt[1] and pt[vi], so points 0 and 1 are on same line
			tempsi = ptorder[planeX].pt[vi];
			ptorder[planeX].pt[vi] = ptorder[planeX].pt[1];
			ptorder[planeX].pt[1] = tempsi;
            
			// calculate cosine between points indexed 1,0,X
			for(vi=2; vi<surfcount; ++vi) {
				cos10X[vi] = cosPQR(poly[ptorder[planeX].pt[1]].xi, poly[ptorder[planeX].pt[0]].xi, 
                                    poly[ptorder[planeX].pt[vi]].xi);
			}
            
			// order by cosines, decreasing order
			for(vi=2; vi<surfcount-1; ++vi) {
				for(vi2=vi+1; vi2<surfcount; ++vi2) {
					if(cos10X[vi] < cos10X[vi2]) {
						// swap indices if points in wrong order
						tempsi = ptorder[planeX].pt[vi];
						ptorder[planeX].pt[vi] = ptorder[planeX].pt[vi2];
						ptorder[planeX].pt[vi2] = tempsi;
						tempcos = cos10X[vi];
						cos10X[vi] = cos10X[vi2];
						cos10X[vi2] = tempcos;
					}
				}
			}
		}
        
	}
	return(surfatom);
}


/**************************
 * subroutine calc_areas
 **************************/

void calc_areas(vertex poly[], const vertex* centerpt, float rado, int NC, int NV, 
                plane cont[], const ptindex* ptorder, const atomsas* atomzero_ptr)
{
	char   engflag;         // ='Y' if an engulfing plane is present
	int    planeX;          // current plane, 1 to NC
	int    NP;              // number of points on face
	int    vi;              // vertices counter
	double area;            // area of polygon
	int    pa1, pa2;        // possible common planes
	int    commplane;       // plane shared by adjacent points (not planeX)
	int    epi;             // engulfing plane counter
	int    engplane[4];     // index to engulfed planes
	double maxSAS;          // maximum solvent exposed surface, no atoms other than engulfing.
	vertex ptB, ptC; // arc point intersections
	int    currpt, nextpt;
	double cosNN1[MAX_CONT];      // angle between vertex N and vertex N+1 
	double cosNzero[MAX_CONT];    // angle between vertex N and vertex N+1 
	double tanprod;         // product of tangents
	int    v0, va, vb;      // vertices for arc calculation
    
	double U,V,W,X;
	double tansqrS, tansqrSA, tansqrSB, tansqrSC;
    
	engflag = 'N';
	epi = 0;
    
	//printf("Calculating AREA for atomzero=[%d]\n",atomzero_ptr->atomnum);
    
	//RESET AREAS TO ZERO
	for(planeX=0; planeX < NC; ++planeX) {
		cont[planeX].area = 0.0;
		if(cont[planeX].flag == 'E') {
			engflag = 'Y';
			engplane[epi] = planeX;
			++epi;
			if (epi==3)
				break;
		}
	}
    
	if(engflag == 'Y') { // engulfing plane correction - project points onto sphere surface.
		project_points(poly, centerpt, rado, NC, NV, cont, atomzero_ptr);
	}
    
	/* ---------------------------- */
	/* calculate area for each face */
	/* ---------------------------- */
    
	for(planeX=0; planeX < NC; ++planeX) {
		NP = ptorder[planeX].numpts;
		area = 0.0;
		if(cont[planeX].flag == 'X') {
			continue;
		}
        
		// if there are no points on a valid contact, area is spherical cap
		if(NP == 0) {
			if(test_point(centerpt[planeX].xi, cont, NC, rado, planeX, -1, -1) == 'Y') {
				cont[planeX].area = 2.0f*PI*rado*(rado-centerpt[planeX].dist);
			}
		} else if(NP == 2) {  // only two contact points, check which part of arc
			if(test_point(centerpt[planeX].xi, cont, NC, rado, planeX, -1, -1) == 'Y') { // area is (cap - arc)
				cont[planeX].area = 2.0f*PI*rado*(rado-centerpt[planeX].dist)
                - spherical_arc(&centerpt[planeX], &poly[ptorder[planeX].pt[0]], 
                                &poly[ptorder[planeX].pt[1]], rado);
			} else {            // area is arc.
				cont[planeX].area = spherical_arc(&centerpt[planeX], &poly[ptorder[planeX].pt[0]], 
                                                  &poly[ptorder[planeX].pt[1]], rado);
			}
		} else {
            
			// ------ calculate cosines and angles ------
			for(vi=0; vi<NP; ++vi) {
				v0 = ptorder[planeX].pt[0];
				va = ptorder[planeX].pt[vi];
				vb = ptorder[planeX].pt[(vi+1)%NP];
                
				// calculate cosines between adjacent vertices
				cosNN1[vi] = (( poly[va].xi[0]*poly[vb].xi[0] + poly[va].xi[1]*poly[vb].xi[1] 
                               + poly[va].xi[2]*poly[vb].xi[2] ) / (poly[va].dist*poly[vb].dist));
                
				// calculate cosines between vertex zero and vertex 'vi'
				if(vi != 0) {
					cosNzero[vi] = ((poly[v0].xi[0]*poly[va].xi[0] + poly[v0].xi[1]*poly[va].xi[1] 
                                     + poly[v0].xi[2]*poly[va].xi[2] ) / (poly[v0].dist*poly[va].dist));
				}
			}
            
			// ----- calculate area of triangles in face -----
			for(vi=1; vi<(NP-1); ++vi) {
				U = sqrt((1+cosNzero[vi])*(1+cosNN1[vi])*(1+cosNzero[vi+1])/8.0);
				V = sqrt((1-cosNzero[vi])*(1-cosNN1[vi])*(1+cosNzero[vi+1])/8.0);
				W = sqrt((1-cosNzero[vi])*(1+cosNN1[vi])*(1-cosNzero[vi+1])/8.0);
				X = sqrt((1+cosNzero[vi])*(1-cosNN1[vi])*(1-cosNzero[vi+1])/8.0);
				tansqrS  = (1-U+V+W+X)/(1+U-V-W-X);
				tansqrSA = (1-U-V-W+X)/(1+U+V+W-X);
				tansqrSB = (1-U-V+W-X)/(1+U+V-W+X);
				tansqrSC = (1-U+V-W-X)/(1+U-V+W+X);
				tanprod = sqrt(tansqrS*tansqrSA*tansqrSB*tansqrSC);
				if(tanprod > 0.0) {
					area += 4.0*rado*rado*atan(sqrt(tanprod));
				}
			}
            
			// ----- add area of arc segments  -----
			for(vi=0; vi<NP; ++vi) {
				va = ptorder[planeX].pt[vi];
				vb = ptorder[planeX].pt[(vi+1)%NP];
                
				//check if adjacent points are arc segments
				if((poly[va].plane[2] == -1) && (poly[vb].plane[2] == -1)) {
					// if on same two planes, no arc.
					if((poly[va].plane[0]+poly[va].plane[1]) != (poly[vb].plane[0]+poly[vb].plane[1])) {
						area += spherical_arc(&centerpt[planeX], &poly[va], &poly[vb], rado);
					}
				}
			}
			cont[planeX].area = area;
		}
	}
    
	// --------------------------------------------------------
	//  add correction terms for engulfing planes, if required
	// --------------------------------------------------------
    
	if(engflag == 'Y') {
		for(planeX=0; planeX < NC; ++planeX) {
			if(cont[planeX].flag != 'E') {
				continue;
			}
            
			NP = ptorder[planeX].numpts;
			for(vi=0; vi<NP; ++vi) {
				currpt = ptorder[planeX].pt[vi];
				nextpt = ptorder[planeX].pt[(vi+1)%NP];
                
				// find common second plane, if any.
				if(poly[currpt].plane[0] == planeX) {
					pa1 = poly[currpt].plane[1];
					pa2 = poly[currpt].plane[2];
				} else {
					pa1 = poly[currpt].plane[0];
					if(poly[currpt].plane[1] == planeX) {
						pa2 = poly[currpt].plane[2];
					} else {
						pa2 = poly[currpt].plane[1];
					}
				}
                
				if((pa1 == poly[nextpt].plane[0]) || (pa1 == poly[nextpt].plane[1])
				   || (pa1 == poly[nextpt].plane[2])) {
					commplane = pa1;
				} else if((pa2 == poly[nextpt].plane[0]) || (pa2 == poly[nextpt].plane[1])
                          || (pa2 == poly[nextpt].plane[2])) { 
					commplane = pa2;
				} else {
					continue;
				}
				if((commplane != -1) && (cont[commplane].flag != 'E')) {
					// add correction to commplane area. here centerpt is from engulfing plane.
					cont[commplane].area += spherical_arc(&centerpt[planeX], &poly[currpt], &poly[nextpt], rado);
					if(NP == 2) break;  // otherwise would repeat adding area
				}
			}
		}
        
		// -----------------------------------------------
		// ------ calculate engulfed contact areas -------
		// -----------------------------------------------
        
		// optimizable
        
		if(epi == 1) {
			cont[engplane[0]].area = 2.0f*PI*rado*(rado+cont[engplane[0]].dist);
		}else if(epi == 2) { 
            
			//printf("engplane[0]=%d engplane[1]=%d\n",engplane[0],engplane[1]);
			//printf("cont[engplane[0]] Ai[0]=%.3f Ai[1]=%.3f Ai[2]=%.3f Ai[3]=%.3f\n",cont[engplane[0]].Ai[0],cont[engplane[0]].Ai[1],cont[engplane[0]].Ai[2],cont[engplane[0]].Ai[3]);
			//printf("cont[engplane[1]] Ai[0]=%.3f Ai[1]=%.3f Ai[2]=%.3f Ai[3]=%.3f\n",cont[engplane[1]].Ai[0],cont[engplane[1]].Ai[1],cont[engplane[1]].Ai[2],cont[engplane[1]].Ai[3]);
            
			//printf("cont[engplane[0]] area=%.3f dist=%.3f\n",cont[engplane[0]].area,cont[engplane[0]].dist);
			//printf("cont[engplane[1]] area=%.3f dist=%.3f\n",cont[engplane[1]].area,cont[engplane[1]].dist);
            
			if(solve_2xS(&cont[engplane[0]],&cont[engplane[1]], rado, ptB.xi, ptC.xi)== -1) {
				cont[engplane[0]].area = 2.0f*PI*rado*rado;
				cont[engplane[1]].area = 2.0f*PI*rado*rado;
			}else {
				ptB.dist = rado;
				ptC.dist = rado;
				maxSAS = spherical_arc(&centerpt[engplane[0]], &ptB, &ptC, rado);
				maxSAS += spherical_arc(&centerpt[engplane[1]], &ptB, &ptC, rado);
				cont[engplane[0]].area = 2.0f*PI*rado*rado - 0.5*maxSAS;
				cont[engplane[1]].area = 2.0f*PI*rado*rado - 0.5*maxSAS;
			}
		} else if(epi>=3) {
			// no exposed surface if there are three or more engulfing contacts 
			for(planeX=0; planeX<NC; ++planeX) {
				if(cont[planeX].flag == 'E') {
					cont[planeX].area = 4.0f*PI*rado*rado/epi;
				} else {
					cont[planeX].area = 0.0f;
				}
			}
		}
	}
	return;
}


/***************************
 * subroutine min_areas
 ***************************/

void min_areas(ca_struct* ca_rec, const atomsas* Calc, const atomsas* atomzero_ptr, char* vcontacts_self_consistency)
{
	int found;
	int currindex, currindex2;
	
	currindex = atomzero_ptr->ca_index;
    
	// loop through contacts of atom A
	while(currindex != -1) {
		
		// for a given contact A --> B
		// if contacts of B were already calculated
		currindex2 = Calc[ca_rec[currindex].atom].ca_index;
		found = 0;
		
		if(currindex2 != -1){
			
			// loops through contacts of atom B
			while(currindex2 != -1){
                
				// B --> A found
				if(ca_rec[currindex].from == ca_rec[currindex2].atom){
                    
					// compare contact areas between A and B
					// FA->vcontacts_self_consistency determines if the mean, min or max contact area is chosen
                    
                    if(!strcmp(vcontacts_self_consistency,"MEAN")){
                        
                        ca_rec[currindex].area = (ca_rec[currindex].area + ca_rec[currindex2].area) / 2.0;
                        ca_rec[currindex2].area = ca_rec[currindex].area;
                        
                        found = 1;
                        break;
                        
                    }else if(!strcmp(vcontacts_self_consistency,"MIN")){
                        
                        if(ca_rec[currindex].area > ca_rec[currindex2].area){
                            ca_rec[currindex].area = ca_rec[currindex2].area;
                        }else{
                            ca_rec[currindex2].area = ca_rec[currindex].area;
                        }
                        
                        found = 1;
                        break;
                        
                    }else if(!strcmp(vcontacts_self_consistency,"MAX")){
                        
                        if(ca_rec[currindex].area < ca_rec[currindex2].area){
                            /*
                             printf("%d --> %d: %.3f now set to %.3f\n", 
                             Calc[ca_rec[currindex].atom].atomnum,
                             Calc[ca_rec[currindex2].atom].atomnum,
                             ca_rec[currindex].area, ca_rec[currindex2].area);
                             */
                            ca_rec[currindex].area = ca_rec[currindex2].area;
                        }else{
                            /*
                             printf("%d --> %d: %.3f now set to %.3f\n", 
                             Calc[ca_rec[currindex2].atom].atomnum,
                             Calc[ca_rec[currindex].atom].atomnum,
                             ca_rec[currindex2].area, ca_rec[currindex].area);
                             */
                            ca_rec[currindex2].area = ca_rec[currindex].area;
                        }
                        
                        found = 1;
                        break;
                        
                    }
                    
                }
                
				currindex2 = ca_rec[currindex2].prev;
			}
            
			// B --> A not found
			// A --> B should then be 0.0
			if(!found){ ca_rec[currindex].area = 0.0f; }
		}
        
		// next contact
		currindex = ca_rec[currindex].prev;
	}
    
}

/***************************
 * subroutine save_areas
 ***************************/

void save_areas(const plane* cont, const contactlist* contlist, int NC, int atomzero,  
                atomsas* Calc, int* ca_recsize,int* numcarec,ca_struct** ca_rec, int* ca_index)
{
	int cai;
	int currindex, previndex, nextindex;
    
	if((*numcarec) > ((*ca_recsize)-100)) {
		(*ca_recsize) += 10000;
		(*ca_rec) = (ca_struct*)realloc((*ca_rec), (*ca_recsize)*sizeof(ca_struct));
		if(!(*ca_rec)) { 
			fprintf(stderr,"ERROR: memory allocation error (*ca_rec)\n"); 
			Terminate(2);
		}
	}
    
	// first overwrite previous records
	cai = 0;
	currindex = ca_index[atomzero];
	previndex = -1;
	while((currindex != -1) && (cai<NC)) {
		if((cont[cai].area > 0.0) && (cont[cai].flag != 'X'))  {
			(*ca_rec)[currindex].from = atomzero;
			(*ca_rec)[currindex].atom = cont[cai].index;
			(*ca_rec)[currindex].area = cont[cai].area;
			(*ca_rec)[currindex].dist = contlist[cai].dist;
			nextindex = (*ca_rec)[currindex].prev; // next index is prev record from old list
			(*ca_rec)[currindex].prev = previndex;
			Calc[atomzero].ca_index = currindex;
			previndex = currindex;
			currindex = nextindex;
		}
		++cai;
	}
	ca_index[atomzero] = previndex;
    
	// then add new records
	while(cai<NC) {
		if((cont[cai].area > 0.0) && (cont[cai].flag != 'X'))  {
			//printf("add new record to atom %d\n", Calc[cont[cai].index].atomnum);
			(*ca_rec)[(*numcarec)].from = atomzero;
			(*ca_rec)[(*numcarec)].atom = cont[cai].index;
			(*ca_rec)[(*numcarec)].area = cont[cai].area;
			(*ca_rec)[(*numcarec)].dist = contlist[cai].dist;
			(*ca_rec)[(*numcarec)].prev = ca_index[atomzero];
			ca_index[atomzero] = (*numcarec);
			Calc[atomzero].ca_index = (*numcarec);
			++(*numcarec);
		}
		++cai;
	}
    
	//printf("atom %d marked as done from save_areas\n", Calc[atomzero].atomnum);
	Calc[atomzero].done = 'Y';
    
	// index contacts for atoms not yet done
	for(cai=0; cai<NC; ++cai) {
		if(Calc[cont[cai].index].done != 'Y') {
			if((cont[cai].area > 0.0) && (cont[cai].flag != 'X'))  {
				//printf("add new record for atom %d (not done yet)\n", Calc[cont[cai].index].atomnum);
				(*ca_rec)[(*numcarec)].from = cont[cai].index;
				(*ca_rec)[(*numcarec)].atom = atomzero;
				(*ca_rec)[(*numcarec)].prev = ca_index[cont[cai].index];
				ca_index[cont[cai].index] = (*numcarec);
				++(*numcarec);
			}
		}
	}
    
	//print_areas(Calc, *numcarec, *ca_rec);
    
	return;
}

/***************************
 * subroutine print_areas
 ***************************/

void print_areas(atomsas* Calc, int numcarec,ca_struct* ca_rec)
{
    
	int i,j;
	int nextatom=-1, atom=-1;
    
	int* done;
	int  ndone=0;
	
	done = (int*)malloc(100000*sizeof(int));
	if(!done){
		printf("ERROR: could not allocate memory for done (print_areas)\n");
		Terminate(1);
	}
    
    while(nextatom){
		atom = nextatom;
		nextatom = -1;
		
		if(atom != -1){ printf("Areas[%5d]: ", Calc[atom].atomnum); }
        
		for(i=0; i<numcarec; i++){
			if(ca_rec[i].from == atom){
				printf("%6d", Calc[ca_rec[i].atom].atomnum);
			}else if(nextatom == -1){
				int f=0;
				for(j=0; j<ndone; j++){
					if(done[j] == ca_rec[i].from){
						f=1;
						break;
					}
				}
				if(!f) { nextatom = ca_rec[i].from; }
			}
		}
		if(atom != -1){	printf("\n"); done[ndone++] = atom; }
		if(nextatom == -1) { nextatom = 0; }
	} 
    
	free(done);
    
	return;
}

/*****************************
 * subroutine project_points
 *****************************/

// this subroutine corrects for engulfing atoms. The center of the engulfing plane
// contact is used as the center of projection instead of the center of the atom.
// points already on the surface are not moved, preserving the SAS.

void project_points(vertex poly[], const vertex* centerpt, float rado, int NC, 
                    int NV, const plane* cont, const atomsas* atomzero_ptr)
{
	int epi;               // engulfing plane counter
	int cai;               // contact atom counter
	int engplane[4];       // index to engulfing planes
	double projpt[3];      // center point for projecting internal pts onto surface of sphere
	double pt0[3], pt1[3]; // temporary points for intersection solutions
	double V[3];           // vector from projection point to vertex
	double *P;             // pointer for vertex coordinates to be projected
	double a,b,c,k;        // for solving quadratic eqn for point projection
	int vi;                // vertex counter
    
	// count and mark engulfing planes
	epi = 0;
	for(cai=0; cai<NC; ++cai) {
		if(cont[cai].flag == 'E') {
			engplane[epi] = cai;
			++epi;
			if(epi==3) 
				break;
		}
	}
    
	// get projpt[] for projecting points to surface
	if(epi == 1) {
		projpt[0] = centerpt[engplane[0]].xi[0];
		projpt[1] = centerpt[engplane[0]].xi[1];
		projpt[2] = centerpt[engplane[0]].xi[2];
	} else if(epi == 2) {
		if (solve_2xS(&cont[engplane[0]],&cont[engplane[1]],rado, pt0, pt1)== -1){
			// added by Francis Gaudreault
			projpt[0] = atomzero_ptr->coor[0] + rado;
			projpt[1] = atomzero_ptr->coor[1];      
			projpt[2] = atomzero_ptr->coor[2];
		} else {
			projpt[0] = (pt0[0]+pt1[0])/2;
			projpt[1] = (pt0[1]+pt1[1])/2;
			projpt[2] = (pt0[2]+pt1[2])/2;
		}
	} else {
		solve_3x3(cont[engplane[0]].Ai, cont[engplane[1]].Ai, cont[engplane[2]].Ai, pt0);
		projpt[0] = pt0[0];
		projpt[1] = pt0[1];
		projpt[2] = pt0[2];
	}
    
	for(vi=0; vi<NV; ++vi) {
		if(poly[vi].plane[2] != -1) {
			// project point to intersection of engulfing plane(s) and surface of sphere
			P = poly[vi].xi;
			V[0] = P[0] - projpt[0];
			V[1] = P[1] - projpt[1];
			V[2] = P[2] - projpt[2];
			a = V[0]*V[0] + V[1]*V[1] + V[2]*V[2];
			b = 2*(P[0]*V[0] +P[1]*V[1] +P[2]*V[2]);
			c = P[0]*P[0] + P[1]*P[1] + P[2]*P[2] - rado*rado;  // c is < 0 
			k = (sqrt(b*b - 4.0*a*c) - b)/(2*a);                // k is > 0
			P[0] += k*V[0];
			P[1] += k*V[1];
			P[2] += k*V[2];      
			poly[vi].dist = rado;
		}
	}
    
	return;
}


/************************
 * function test_point     
 ************************/

// this function tests a given point versus a set of plane constraints
// it returns a value of 'Y' if the point is OK, and 'N' if there
// was a violation of the plane inequality.
// ptX =(xo, yo, zo); if Axo+Byo+Czo+D > 0, point is behind plane (hidden)
// planes A,B,C are the planes that the point lies on, don't test.

char test_point(const double* ptX, const plane* cont, int NC, float rado, int planeA, int planeB, int planeC)
{
	int cp;      // counter, number of planes
    
	// if point is not behind any plane, keep point.
	for(cp=0; cp<NC; ++cp) {
		//if pt is on current plane, get next plane
		if((cp != planeA) && (cp != planeB) && (cp != planeC) && (cont[cp].flag != 'X')) {
			if((cont[cp].Ai[0]*ptX[0] + cont[cp].Ai[1]*ptX[1] + cont[cp].Ai[2]*ptX[2] + cont[cp].Ai[3]) > 0.0) {
				//point is behind plane, not on polyhedron.
				return('N');
			}
		}
	}
	return('Y');
}


/****************************
 * function solve_3x3   
 ****************************/

// determines the intersection of three planes
// (solves a system of three linear equations and 3 unknowns)
// input is 3 four element arrays, representing Ax+By+Cz+D=0
// output is a three element array, (xi,xj,xk).

int solve_3x3(const double* eq0, const double* eq1, const double* eq2, double pt[])
{
	double cof00, cof01, cof02;   // matrix cofactors
	double cof10, cof11, cof12;
	double cof20, cof21, cof22;
	double det;                   // determinant of matrix
    
	cof00 =  eq1[1]*eq2[2] - eq2[1]*eq1[2];
	cof01 = -eq1[0]*eq2[2] + eq2[0]*eq1[2];
	cof02 =  eq1[0]*eq2[1] - eq2[0]*eq1[1];
    
	cof10 = -eq0[1]*eq2[2] + eq2[1]*eq0[2];
	cof11 =  eq0[0]*eq2[2] - eq2[0]*eq0[2];
	cof12 = -eq0[0]*eq2[1] + eq2[0]*eq0[1];
    
	cof20 =  eq0[1]*eq1[2] - eq1[1]*eq0[2];
	cof21 = -eq0[0]*eq1[2] + eq1[0]*eq0[2];
	cof22 =  eq0[0]*eq1[1] - eq1[0]*eq0[1];
    
	det = eq0[0]*cof00 + eq0[1]*cof01 + eq0[2]*cof02;
	if(det == 0.0) {
		//printf("no solution for equation set, determinant is zero\n");
		return(-1);
	} else {
		pt[0] = -(eq0[3]*cof00 + eq1[3]*cof10 + eq2[3]*cof20)/det;
		pt[1] = -(eq0[3]*cof01 + eq1[3]*cof11 + eq2[3]*cof21)/det;
		pt[2] = -(eq0[3]*cof02 + eq1[3]*cof12 + eq2[3]*cof22)/det;
		//printf("3x3 pt: 0=[%5.3f]\t1=[%5.3f]\t2=[%5.3f]\n",pt[0],pt[1],pt[2]);
		return(0);
	}
}

/**********************************
 * Function solve_2xS       
 * revised 19/02/2001  BJM
 ***********************************/

// determines the intersection of two planes and a sphere radius 'rad'
// input is 2 four element arrays, representing Ax+By+Cz+D=0
// plus the radius of a sphere centered on the origin.
// output is two three element arrays pt0 and pt1, (xi,xj,xk).
// return value is -1 if no real solution exists.

int solve_2xS(const plane* eq0,const plane* eq1, float rado, double pt0[], double pt1[])
{
	double eq2[3];              // eqn of plane through (0,0,0) and perp. to other two
	double cof00, cof01, cof02; // matrix cofactors
	double cof10, cof11, cof12; // (don't need cof20, cof21, cof22)
	double det;                 // determinant of matrix
	double avgpt[3];            // average of two solution points
	double t;                   // parameter in eqn of line: x=xo+At, y=yo+Bt, z=zo+Ct. 
	int xi;                     // coordinate counter
    
	eq2[0] = eq0->Ai[1]*eq1->Ai[2] - eq0->Ai[2]*eq1->Ai[1];
	eq2[1] = eq0->Ai[2]*eq1->Ai[0] - eq0->Ai[0]*eq1->Ai[2];
	eq2[2] = eq0->Ai[0]*eq1->Ai[1] - eq0->Ai[1]*eq1->Ai[0];
    
	cof00 =  eq1->Ai[1]*eq2[2] - eq2[1]*eq1->Ai[2];
	cof01 = -eq1->Ai[0]*eq2[2] + eq2[0]*eq1->Ai[2];
	cof02 =  eq1->Ai[0]*eq2[1] - eq2[0]*eq1->Ai[1];
    
	cof10 = -eq0->Ai[1]*eq2[2] + eq2[1]*eq0->Ai[2];
	cof11 =  eq0->Ai[0]*eq2[2] - eq2[0]*eq0->Ai[2];
	cof12 = -eq0->Ai[0]*eq2[1] + eq2[0]*eq0->Ai[1];
    
	det = eq2[0]*eq2[0] + eq2[1]*eq2[1] + eq2[2]*eq2[2];
	if(det == 0.0) {
		//printf("no solution in solve_2xS\n");
		return(-1);
	}
    
	avgpt[0] = -(eq0->Ai[3]*cof00 + eq1->Ai[3]*cof10)/det;
	avgpt[1] = -(eq0->Ai[3]*cof01 + eq1->Ai[3]*cof11)/det;
	avgpt[2] = -(eq0->Ai[3]*cof02 + eq1->Ai[3]*cof12)/det;
    
	t = (rado*rado-avgpt[0]*avgpt[0]-avgpt[1]*avgpt[1]-avgpt[2]*avgpt[2])/det;
	if(t<0.0) {
		//printf("t smaller than 0\n");
		return(-1);
	} else {
		t = sqrt(t);
	}
    
	for(xi=0; xi<3; ++xi) {
		pt0[xi] = avgpt[xi] + t*eq2[xi];
		pt1[xi] = avgpt[xi] - t*eq2[xi];
	}
    
	//printf("ptB.xi 0=[%.3f] 1=[%.3f] 2=[%.3f]\n",pt0[0],pt0[1],pt0[2]);
	//printf("ptC.xi 0=[%.3f] 1=[%.3f] 2=[%.3f]\n",pt1[0],pt1[1],pt1[2]);
    
	return(0);
}


/***************************
 * function cosPQR         *
 ***************************/

// this function returns the cosine of the angle between three points P,Q,R 
// with Q being the center point.

double cosPQR(const double* ptP,const double* ptQ,const double* ptR)
{
	double QP[3];    // vector from Q to P
	double QR[3];    // vector from Q to R
	double cosine;   // cosine of angle PQR at Q.
    
	// calculate vectors
	QP[0] = ptP[0] - ptQ[0];
	QP[1] = ptP[1] - ptQ[1];
	QP[2] = ptP[2] - ptQ[2];
	QR[0] = ptR[0] - ptQ[0]; 
	QR[1] = ptR[1] - ptQ[1]; 
	QR[2] = ptR[2] - ptQ[2]; 
    
	//calculate cosine
	cosine = (QP[0]*QR[0]+ QP[1]*QR[1]+ QP[2]*QR[2])
    /sqrt((QP[0]*QP[0]+ QP[1]*QP[1]+ QP[2]*QP[2]) * (QR[0]*QR[0]+ QR[1]*QR[1]+ QR[2]*QR[2]));
    
	return(cosine);
}


/********************************
 * function spherical_arc       *	   
 ********************************/

// given two points in cartesian space and the center of the spherical
// cap between atom A and the origin, plus the radius of the sphere,
// this function returns the area of an arc between points B and C.
// the sides of the arc are the great circle distance (shortest distance)
// and the arc of the spherical cap centered on line AO.

double spherical_arc(const vertex* ptAo,const vertex* ptB,const vertex* ptC, float rado)
{
	// here, cosAOC=cosAOB, AOB = AOC.
	// BAC is the angle at the point on the origin, on line AO
    
	double cosAOB, cosBOC; // cosines of angles centered on (0,0,0)
	double cosBAC, angBAC;    // angle and cosine at vertex of cap
	double U, V, W;
	double tansqrS, tansqrSA, tansqrSB;
	double tanprod;        // product of tangents for spherical triangle
	double area;           // the value to be returned
    
	cosAOB = (ptAo->xi[0]*ptB->xi[0] + ptAo->xi[1]*ptB->xi[1] + ptAo->xi[2]*ptB->xi[2])/(ptAo->dist*ptB->dist);
	cosBOC = (ptB->xi[0]*ptC->xi[0] + ptB->xi[1]*ptC->xi[1] + ptB->xi[2]*ptC->xi[2])/(ptB->dist*ptC->dist);
    
	U = (1+cosAOB)*sqrt((1+cosBOC)/8.0);
	V = (1-cosAOB)*sqrt((1+cosBOC)/8.0);
	W = sqrt((1-cosAOB)*(1+cosAOB)*(1-cosBOC)/8.0);  // W == X
	tansqrS  = (1-U+V+W+W)/(1+U-V-W-W);
	tansqrSB = (1-U-V)/(1+U+V);
	tansqrSA = (1-U+V-W-W)/(1+U-V+W+W);
	if((tansqrS*tansqrSA) > 0.0) {
		tanprod = sqrt(tansqrSB*tansqrSB*tansqrS*tansqrSA);
	} else {
		tanprod = 0.0f;
	}
    
	cosBAC = cosPQR(ptB->xi, ptAo->xi, ptC->xi);
	if(cosBAC>1.0) { 
		angBAC = 0.0f; 
	} else if(cosBAC<-1.0) {
		angBAC = PI; 
	} else { 
		angBAC = acos(cosBAC);
	} 
    
	// area is area of spherical cap segment minus area of spherical triangle
	area = rado*(angBAC*(rado-ptAo->dist) - 4.0*rado*atan(sqrt(tanprod)));
	return(area);
}

/******************************
 * subroutine index_protein
 * created 16/02/2001
 ******************************
 
 // assigns all protein atoms to boxes within a cubic grid
 // returns cube dimesions as number of boxes per side
 // assigns values to global array box[]
 
 *****************************/

void index_protein(FA_Global* FA,atom* atoms,resid* residue,atomsas* Calc,atomindex** box,int* Calclist,int* dim,int atmcnt)
{
	int   i,j;             // dumb counters
	int   resi;            // residue counter
	int   atmi;            // atom counter 0 to (atmcnt-1)
	int   boxi;            // box index
	int   startind;        // start point of box in Calclist
	int   rot;             // sc rotamer
	float max_width;
	float diff;
	int   alter;
	float global_min[3];
	float global_max[3];
	int   dim2,dim3;
    
	rot=0;
	i=0;
	alter=0;
    
	// Do not alter default PDB min and max coordinates
	for(j=0;j<3;++j){
		global_min[j]=FA->globalmin[j];
		global_max[j]=FA->globalmax[j];
	}
    
	max_width=FA->maxwidth;
    
    
	for (resi=1; resi<=FA->res_cnt; ++resi) {
		rot = residue[resi].rot;
        
		//printf("-----residue[%d]=%s - rot=%d-----\n",residue[resi].number,residue[resi].name,residue[resi].rot);
        
		for (atmi=residue[resi].fatm[rot];atmi<=residue[resi].latm[rot];++atmi){
			// Copy atoms structure to the new Vcont structure
			// only the atoms that correspond to the correct rotamer are copied
			// the total number of atoms thus is equal to atm_cnt_real
            
			Calc[i].atomnum=atoms[atmi].number;
			for (j=0;j<3;j++){
				Calc[i].coor[j]=atoms[atmi].coor[j];
			}
            
			strcpy(Calc[i].atomname,atoms[atmi].name);
			strcpy(Calc[i].res,residue[resi].name);
			Calc[i].resnum=residue[resi].number;
			Calc[i].inum=atoms[atmi].ofres;
			Calc[i].chn=residue[resi].chn;
			Calc[i].radius=atoms[atmi].radius;
			Calc[i].type=atoms[atmi].type;
			Calc[i].done='N';
			Calc[i].ca_index=-1;
			Calc[i].score=(atoms[atmi].optres != NULL);
            
            
			for(j=0;j<3;++j){
				if (atoms[atmi].coor[j] < global_min[j]){
					global_min[j]=atoms[atmi].coor[j];
					alter=1;
				}
				if (atoms[atmi].coor[j] > global_max[j]){
					global_max[j]=atoms[atmi].coor[j];
					alter=1;
				}
			}
            
			++i;
		}
	}
    
	//printf("[%d] atoms were copied to Vcont atomsas_struct compared to real [%d]\n", i, FA->atm_cnt_real);
    
    
	if (alter) {
		for(j=0;j<3;++j){
			diff=global_max[j]-global_min[j];
			if (diff > max_width) {
				//printf("New difference=%8.3f\n",diff);
				max_width=diff;
			}
		}
	}
    
	// ------ get largest dimension of protein -------
	*dim = 0;
	*dim = (int)(max_width/CELLSIZE)+1;
	//printf("New Dimension=[%d]\n",dim);
    
	dim2 = (*dim)*(*dim);
	dim3 = (*dim)*(*dim)*(*dim);
    
	(*box) = (atomindex*)malloc(dim3*sizeof(atomindex));
	if(!(*box)){
		fprintf(stderr,"ERROR: memory allocation error for box\n");
		Terminate(2);
	}
    
	memset((*box),0,dim3*sizeof(atomindex));  
    
	// count entries per box, assign box number to atom
	for(atmi=0;atmi<atmcnt;++atmi){
		boxi = (int)((Calc[atmi].coor[0]-global_min[0])/CELLSIZE)*dim2
        + (int)((Calc[atmi].coor[1]-global_min[1])/CELLSIZE)*(*dim)
        + (int)((Calc[atmi].coor[2]-global_min[2])/CELLSIZE);
		Calc[atmi].boxnum = boxi;
		//printf("coor[0]: %8.3f\tcoor[1]: %8.3f\tcoor[2]: %8.3f\n", Calc[atmi].coor[0],Calc[atmi].coor[1],Calc[atmi].coor[2]);
		//printf("Calc[%d]=%d Boxi=[%d] RNum=[%d]\n",atmi,Calc[atmi].atomnum,boxi,Calc[atmi].resnum);
        
		++(*box)[Calc[atmi].boxnum].nument;
	}
    
    
	// assign start pointers for boxes in Calclist
	startind = 0;
	for (boxi=0; boxi<dim3; ++boxi) {
		(*box)[boxi].first = startind;
		startind += (*box)[boxi].nument;
	}
    
	// clear array (needed for recounting index)
	for(boxi=0; boxi<dim3; ++boxi) {
		(*box)[boxi].nument = 0;
	}
    
	// fill Calclist index
	for (atmi=0; atmi<atmcnt; ++atmi) {
		boxi = Calc[atmi].boxnum;
		Calclist[(*box)[boxi].first+(*box)[boxi].nument] = atmi;
		++(*box)[boxi].nument;
	}
    
}

/*********************************
 * subroutine get_contlist4
 *********************************/

// uses box index of protein atoms to find contacts in range of a given atom.
// requires global variable 'box[]'.
// checks previous atoms, keeps only those with non-zero contact area.

int get_contlist4(atom* atoms,int atomzero, contactlist contlist[], 
                  int atmcnt, float rado, int dim, atomsas* Calc, 
                  const int* Calclist, const atomindex* box, 
                  const ca_struct* ca_rec,const int* ca_index) 
{
	double sqrdist;             // distance squared between two points
	double neardist;            // max distance for contact between atom spheres
	int NC;                     // number of contacts
	int bai;                    // box atom counter
	int boxi;                   // current box number
	int atomj;                  // index number of atom in Calc
	int boxzero;                // box atomzero is in
	int i;                      // dummy counter for surrounding boxes
	int currindex;
	char recalc;                // recalculate neglecting done atoms
    
	int dim2,dim3;
    
	//printf("=====ATOMZERO[%d]=====\n",Calc[atomzero].atomnum);
    
	NC = 0;
	recalc = 'N';
    
	dim2 = dim*dim;
	dim3 = dim*dim*dim;
    
	// mark previously contacted atoms
	currindex = ca_index[atomzero];
	while(currindex != -1) {
		//printf("atom %d marked as contact\n",Calc[ca_rec[currindex].atom].atomnum);
		Calc[ca_rec[currindex].atom].done = 'C'; // makes contact
		currindex = ca_rec[currindex].prev;
	}
    
	// get pdb atom contacts from current and adjacent boxes
	boxzero = Calc[atomzero].boxnum;
	//printf("boxzero=%d\n",boxzero);
	for(i=0; i<27; ++i) {
		// get up to 27 boxes surrounding current box
		boxi = boxzero +dim2*((i/9)-1) + dim*(((i/3)%3)-1) +(i%3)-1;
		if((boxi < 0) || (boxi >= dim3)) continue;  // don't do boxes outside of range
        
		//printf("boxi=%d\n",boxi);
		bai = 0;
		while(bai<box[boxi].nument) {
			//printf("nument=%d\tbai=%d\n",box[boxi].nument,bai);
            
			atomj = Calclist[box[boxi].first+bai]; 
            
			// check previous contacts
			// if(recalc == 'Y' && Calc[atomj].done == 'Y') {
            
			if(Calc[atomj].done == 'Y') {
				//printf("skipped atom %d\n",Calc[atomj].atomnum);
				++bai;
				continue;
			}
            
            
			sqrdist = (Calc[atomzero].coor[0]-Calc[atomj].coor[0])*(Calc[atomzero].coor[0]-Calc[atomj].coor[0]) 
            + (Calc[atomzero].coor[1]-Calc[atomj].coor[1])*(Calc[atomzero].coor[1]-Calc[atomj].coor[1])
            + (Calc[atomzero].coor[2]-Calc[atomj].coor[2])*(Calc[atomzero].coor[2]-Calc[atomj].coor[2]);
			neardist =  rado + Calc[atomj].radius + Rw;
            
            
			//printf("neardist = rado(%.3f) + atomj.rad(%.3f) + Rw(%.3f)\n",rado,Calc[atomj].radius,Rw);
			//printf("atom %d is sqrdist(%5.2fA) & neardist(%5.2fA) from atom %d\n",Calc[atomzero].atomnum,sqrdist,neardist*neardist,Calc[atomj].atomnum);
            
            
			if((sqrdist < neardist*neardist) && (sqrdist != 0.0)) {
                
				// add atoms to list
				//printf("atom %d is in contact with atom %d (%.3f)...\n",Calc[atomzero].atomnum,Calc[atomj].atomnum,sqrtf(sqrdist));
                
				contlist[NC].index = atomj;
				contlist[NC].dist = sqrt(sqrdist);
				++NC;
			}
			++bai;
		}
	}
    
	/*
     if (NC==0) {
     recalc = 'Y';
     FA->recalci++;
     //PAUSE;
     goto RESTART;
     }
     */
    
	// reset atoms to 'done'
	currindex = ca_index[atomzero];
	while(currindex != -1) {
		//printf("atom %d marked as done\n",Calc[ca_rec[currindex].atom].atomnum);
		Calc[ca_rec[currindex].atom].done = 'Y'; // re-mark as done
		currindex = ca_rec[currindex].prev;
	}
    
	return(NC);
    
}
