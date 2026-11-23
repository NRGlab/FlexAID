#include "flexaid.h"
#include "boinc.h"

void assign_eigen(FA_Global *FA,atom* atoms,resid* residue, int res_cnt, int nmodes){
    int i,j,k,l;
    
    for(i=1;i<=res_cnt;i++){
        
        // Do not apply eigenvectors on HET groups
        if(residue[i].type == 1)
            continue;
        
        for(j=residue[i].fatm[0];j<=residue[i].latm[0];j++){
            
            atoms[j].eigen = (float**)malloc(nmodes*sizeof(float*));
            if (!atoms[j].eigen){
                fprintf(stderr,"ERROR: memory allocation error for eigen\n");
                Terminate(2);
            }      
            memset(atoms[j].eigen,NULL,nmodes*sizeof(float*));
            
            for(k=0;k<nmodes;k++){
                atoms[j].eigen[k] = (float*)malloc(3*sizeof(float));
                if(!atoms[j].eigen[k]){
                    fprintf(stderr,"ERROR: memory allocation error for eigen\n");
                    Terminate(2);
                }
                memset(atoms[j].eigen[k],0,3*sizeof(float));
            }
            
            for(k=0;k<nmodes;k++){
                
                if(FA->supernode){
                    // N or H
                    if ((strncmp(atoms[j].name," N  ",4) == 0) || (strncmp(atoms[j].name," H  ",4) == 0)){
                        for(l=0;l<3;l++){
                            atoms[j].eigen[k][l]=FA->eigenvector[(i-1)*9+l][k];
                            //printf("assigned %.4f to atom[%d]\n",atoms[j].eigen[k][l],atoms[j].number);
                        }
                        // C, O or OXT
                    }else if ((strncmp(atoms[j].name," C  ",4) == 0) || (strncmp(atoms[j].name," O  ",4) == 0) ||
                              (strncmp(atoms[j].name," OXT",4) == 0)){
                        for(l=0;l<3;l++){
                            atoms[j].eigen[k][l]=FA->eigenvector[(i-1)*9+6+l][k];
                            //printf("assigned %.4f to atom[%d]\n",atoms[j].eigen[k][l],atoms[j].number);
                        }
                    }else{ // side-chain atoms + CA + HA
                        for(l=0;l<3;l++){
                            atoms[j].eigen[k][l]=FA->eigenvector[(i-1)*9+3+l][k];
                            //printf("assigned %.4f to atom[%d]\n",atoms[j].eigen[k][l],atoms[j].number);
                        }
                    }
                }else{
                    
                    // all atoms from a residue have the same eigenvectors
                    for(l=0;l<3;l++){
                        atoms[j].eigen[k][l]=FA->eigenvector[(i-1)*3+l][k];
                        //printf("assigned %.4f to atom[%d]\n",atoms[j].eigen[k][l],atoms[j].number);
                    }
                }
            }
        }
    }
    
    return;
}
