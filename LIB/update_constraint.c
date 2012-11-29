#include "flexaid.h"
#include "boinc.h"

void update_constraint(atom* atoms, int index, constraint* cons)
{
  
	if(atoms[index].ncons == 0){
		atoms[index].cons = (constraint**)malloc(sizeof(constraint*));
	}else{
		atoms[index].cons = (constraint**)realloc(atoms[index].cons,(atoms[index].ncons+1)*sizeof(constraint*));
	}
	
	if(atoms[index].cons == NULL){
		fprintf(stderr,"ERROR: Could not allocate memory for atoms[%d].cons\n",index);
		Terminate(2);
	}
	
	atoms[index].cons[atoms[index].ncons] = cons;
	
	/*
	  printf("atoms[%d].cons[%d] points to cons[%p].\n",
	  index,atoms[index].ncons,cons);
	*/
	
	/*
	  printf("new constraint for atom[%d] with cons.id[%d].\n",
	  atoms[index].number,atoms[index].cons[atoms[index].ncons]->id);
	*/
	
	++atoms[index].ncons;
	
}
