#include "flexaid.h"

int dee_pivot(psFlexDEE_Node psFlexDEEInsertNode,psFlexDEE_Node* psFlexDEENode, int sta, int end, int pivot, int pos, int num_rot)
{
  
  int i,j;
  
  int new_sta = -1;
  int new_end = -1;
  int new_pivot;
  
  //printf("dee_pivot[%d]\tsta=%d\tend=%d\tpivot=%d\tpos=%d  ",x++,sta,end,pivot,pos);

  //if( x % 10 == 0 ) { getchar(); }

  i = pivot - pos;

  if( i < 0 ){
    
    while( i != 0 ){
      (*psFlexDEENode) = (*psFlexDEENode)->prev;
      i++;
    }
    
  }else if( i > 0 ){
    
    while( i != 0 ){
      (*psFlexDEENode) = (*psFlexDEENode)->next;
      i--;
    }
    
  }

  j = cmp_FlexDEE_Nodes(psFlexDEEInsertNode,(*psFlexDEENode),num_rot);
  
  //printf("   ==> %s\n", j < 0?"OVER":"BELOW");

  if ( j == 0 ) { return j; }

  if ( (sta - end) >= 0 ) {
    //printf("   END ==> %s\n", j > 0?"OVER":"BELOW");
    
    return j;
  }

  
  if ( j == -1 ) {
    new_sta = sta;
    new_end = pivot - 1;
  }else if ( j == 1 ) {
    new_sta = pivot + 1;
    new_end = end;
  }
  
  new_pivot = (int)(new_sta + (new_end - new_sta)/2);

  // impossible to have the same rotlist twice
  //if( j == 0 ){ return psFlexDEENode->prev; }

  return dee_pivot(psFlexDEEInsertNode,psFlexDEENode,new_sta,new_end,new_pivot,pivot,num_rot);

}


int cmp_FlexDEE_Nodes(psFlexDEE_Node Node1, psFlexDEE_Node Node2, int num_rot)
{

  int i;

  for(i=0;i<num_rot;++i){

    if(Node1->rotlist[i] > Node2->rotlist[i]){
      return 1;
    }else if(Node1->rotlist[i] < Node2->rotlist[i]){
      return -1;
    }else{
      continue;
    }
    
  }
  
  return 0;
}


void dee_first(psFlexDEE_Node Node, psFlexDEE_Node FirstNode)
{

  while ( Node->prev != NULL ) {

    Node = Node->prev;
    
  }

  while ( Node->next != NULL ) {

    Node->first = FirstNode;
    Node = Node->next;
    
  }
  
  Node->first = FirstNode;
  
  return;

}



void dee_last(psFlexDEE_Node Node, psFlexDEE_Node LastNode)
{

  while ( Node->prev != NULL ) {

    Node = Node->prev;
    
  }

  while ( Node->next != NULL ) {

    Node->last = LastNode;
    Node = Node->next;
    
  }
  
  Node->last = LastNode;
  
  return;

}
