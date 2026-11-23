#include "flexaid.h"

void dee_print(psFlexDEE_Node psFlexDEENode, int num_rot)
{

  int i,j;

  printf("DEE sorted list:\n");
  
  if ( psFlexDEENode == NULL ) { printf("<EMPTY>\n"); return; }


  // loop through first item
  while ( psFlexDEENode->prev != NULL ) {
    
    psFlexDEENode = psFlexDEENode->prev;

  }

  
  j = 0 ;
  
  do {
    
    printf("[%d]: ",++j);
    
    for(i=0;i<num_rot;i++){
      printf("%3d",psFlexDEENode->rotlist[i]);
    }

    printf("\n");

    if ( psFlexDEENode->next != NULL ) { 
      psFlexDEENode = psFlexDEENode->next; 
      
    } else { break; }

  } while ( 1 );

  //printf("rotlist[%d].next is NULL = ",j);for(i=0;i<4;i++){printf("%3d",psFlexDEENode->rotlist[i]);}printf("\n");

  return;

}
