#include "flexaid.h"

// This function assigns an internal probability to each flexible side-chain
// You can override these internal prob. by adding extra PROBLT lines in the flexscfile 
// The internal probabilities come from the analysis of the Holo-Apo Protein Pairs (HAP2db)

void set_intprob(flxsc* flex_res)
{
  
  if(!strcmp(flex_res[0].name,"VAL")){
    flex_res[0].prob=0.070f;
  }else if(!strcmp(flex_res[0].name,"LEU")){
    flex_res[0].prob=0.136f;
  }else if(!strcmp(flex_res[0].name,"ILE")){
    flex_res[0].prob=0.172f;
  }else if(!strcmp(flex_res[0].name,"MET")){
    flex_res[0].prob=0.342f;
  }else if(!strcmp(flex_res[0].name,"ASN")){
    flex_res[0].prob=0.279f;
  }else if(!strcmp(flex_res[0].name,"CYS")){
    flex_res[0].prob=0.048f;
  }else if(!strcmp(flex_res[0].name,"SER")){
    flex_res[0].prob=0.195f;
  }else if(!strcmp(flex_res[0].name,"THR")){
    flex_res[0].prob=0.058f;
  }else if(!strcmp(flex_res[0].name,"GLN")){
    flex_res[0].prob=0.347f;
  }else if(!strcmp(flex_res[0].name,"ASP")){
    flex_res[0].prob=0.220f;
  }else if(!strcmp(flex_res[0].name,"GLU")){
    flex_res[0].prob=0.188f;
  }else if(!strcmp(flex_res[0].name,"LYS")){
    flex_res[0].prob=0.532f;
  }else if(!strcmp(flex_res[0].name,"ARG")){
    flex_res[0].prob=0.381f;
  }else if(!strcmp(flex_res[0].name,"HIS")){
    flex_res[0].prob=0.168f;
  }else if(!strcmp(flex_res[0].name,"PHE")){
    flex_res[0].prob=0.203f;
  }else if(!strcmp(flex_res[0].name,"TRP")){
    flex_res[0].prob=0.090f;
  }else if(!strcmp(flex_res[0].name,"TYR")){
    flex_res[0].prob=0.190f;
  }else{
    flex_res[0].prob=0.000f;
  }


  return;
}
