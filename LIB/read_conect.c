#include "flexaid.h"
/******************************************************************************
 * SUBROUTINE read_conect reads the conectivity data present in the PDB file
 * in case such information is present.
 ******************************************************************************/
void read_conect(FA_Global* FA,atom** atoms,char line[]){
    
    //char num_char[6];                 /* string to read integer */
    int i,j;                          /* dumb counter */
    int atm;                          /* number of the atom we read info */
    //int atn;                          /* atn is covelently bonded to atm */
    int tmp[4];
    char number[6];
    unsigned int len;
    int ncont;
    
    //CONECT LINE
    //CONECT%5d%5d\n
    len = (int)(strlen(line));
    //printf("Connection line=%u\n",len);
    len -= 12; //Remove 'CONECT' and '\n' and First Atom '%5d'
    ncont = len/5;
    
    for(i=0;i<=3;i++){tmp[i]=0;}
    
    for(i=0;i<5;i++){
        number[i]=line[i+6];
    }
    number[5]='\0';
    sscanf(number,"%d",&atm);
    
    //printf("Atom number[%d] has %d contact(s)\n",atm,ncont);
    for(i=0;i<ncont;i++){
        for(j=0;j<5;j++){
            number[j]=line[i*5+j+11];
        }
        number[5]='\0';
        sscanf(number,"%d",&tmp[i]);
    }
    
    //printf("CONECT atm=%d tmp[0]=%d tmp[1]=%d tmp[2]=%d tmp[3]=%d\n",atm,tmp[0],tmp[1],tmp[2],tmp[3]);
    
    for(i=0;i<=3;i++){
        if(tmp[i] != 0){
            //printf("Added for atom[%d]:%d\n",atm,tmp[i]);
            (*atoms)[FA->num_atm[atm]].bond[i+1]=FA->num_atm[tmp[i]];
            (*atoms)[FA->num_atm[atm]].bond[0]++;
        }else{
            break;
        }
    }
    
    return;
}
