#include "flexaid.h"

// This function assigns an internal probability to each flexible side-chain
// You can override these internal prob. by adding extra PROBLT lines in the flexscfile 
// The internal probabilities come from the analysis of the Holo-Apo Protein Pairs (HAP2db)

void set_intprob(flxsc* flex_res)
{
  
  if(!strcmp(flex_res[0].name,"VAL")){
    flex_res[0].prob=0.075f;
  }else if(!strcmp(flex_res[0].name,"LEU")){
    flex_res[0].prob=0.105f;
  }else if(!strcmp(flex_res[0].name,"ILE")){
    flex_res[0].prob=0.159f;
  }else if(!strcmp(flex_res[0].name,"MET")){
    flex_res[0].prob=0.293f;
  }else if(!strcmp(flex_res[0].name,"ASN")){
    flex_res[0].prob=0.233f;
  }else if(!strcmp(flex_res[0].name,"CYS")){
    flex_res[0].prob=0.091f;
  }else if(!strcmp(flex_res[0].name,"SER")){
    flex_res[0].prob=0.180f;
  }else if(!strcmp(flex_res[0].name,"THR")){
    flex_res[0].prob=0.077f;
  }else if(!strcmp(flex_res[0].name,"GLN")){
    flex_res[0].prob=0.293f;
  }else if(!strcmp(flex_res[0].name,"ASP")){
    flex_res[0].prob=0.168f;
  }else if(!strcmp(flex_res[0].name,"GLU")){
    flex_res[0].prob=0.193f;
  }else if(!strcmp(flex_res[0].name,"LYS")){
    flex_res[0].prob=0.440f;
  }else if(!strcmp(flex_res[0].name,"ARG")){
    flex_res[0].prob=0.289f;
  }else if(!strcmp(flex_res[0].name,"HIS")){
    flex_res[0].prob=0.129f;
  }else if(!strcmp(flex_res[0].name,"PHE")){
    flex_res[0].prob=0.183f;
  }else if(!strcmp(flex_res[0].name,"TRP")){
    flex_res[0].prob=0.075f;
  }else if(!strcmp(flex_res[0].name,"TYR")){
    flex_res[0].prob=0.165f;
  }else{
    flex_res[0].prob=0.000f;
  }


  return;
}
